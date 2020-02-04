"""Value Reefs."""
import argparse
import glob
import logging
import os
import sys

import ecoshard
from osgeo import gdal
import pygeoprocessing
import taskgraph
import taskgraph_downloader_pnn


ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')

GLOBAL_REEFS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_reef_md5_5a90d55a505813b5aa9662faee351bf8.tif')

LS_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')
POVERTY_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'Poverty_Count_2017_clean_md5_0cc6e0187be07e760e66f759a0a1f7e8.tif')

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

WORKSPACE_DIR = 'value_reefs_workspace'
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')


def align_raster_list(raster_path_list, target_directory):
    """Aligns all the raster paths.

    Rasters are aligned using the pixel size of the first raster and use
    the intersection and near interpolation methods.

    Parameters:
        raster_path_list (list): list of str paths to rasters.
        target_directory (str): path to a directory to hold the aligned
            rasters.

    Returns:
        list of raster paths that are aligned with intersection and near
            interpolation algorithm.

    """
    if not hasattr(align_raster_list, 'task_graph_map'):
        align_raster_list.task_graph_map = {}
    if target_directory not in align_raster_list.task_graph_map:
        align_raster_list.task_graph_map[target_directory] = (
            taskgraph.TaskGraph(target_directory, -1))
    task_graph = align_raster_list.task_graph_map[target_directory]
    aligned_path_list = [
        os.path.join(target_directory, os.path.basename(path))
        for path in raster_path_list]
    target_pixel_size = pygeoprocessing.get_raster_info(
        raster_path_list[0])['pixel_size']
    LOGGER.debug('about to align: %s', str(raster_path_list))
    task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            raster_path_list, aligned_path_list,
            ['near'] * len(raster_path_list), target_pixel_size,
            'intersection'),
        target_path_list=aligned_path_list)
    return aligned_path_list


def calculate_habitat_population_value(
        shore_sample_point_vector_path, population_raster_path_id_list,
        dem_raster_path, habitat_fieldname_list, habitat_vector_path_map,
        results_dir):
    """Calculate population within protective range of habitat.

    Parameters:
        shore_sample_point_vector_path (str): path to a point shapefile that is only
            used for referencing the points of interest on the coastline.
        population_raster_path_id_list (list): list of (raster_path, field_id)
            tuples. The values in the raster paths will be masked where it
            overlaps with < 10m dem height and convolved within 2km. That
            result is in turn spread onto the habitat coverage at a distance
            of the protective distance of that habitat. These rasters are in
            wgs84 lat/lng projection.
        dem_raster_path (str): path to a dem used to mask population by height
            in wgs84 lat/lng projection.
        habitat_fieldname_list (list): list of habitat ids to analyse.
        habitat_vector_path_map (dict): maps fieldnames from
            `habitat_fieldname_list` to 3-tuples of
            (path to hab raster (str), risk val (float),
             protective distance (float)).
        results_dir (str): path to directory containing habitat back projection
            results

    Returns:
        None

    """
    temp_workspace_dir = os.path.join(
        results_dir, 'calc_pop_coverage_churn')
    taskgraph_working_dir = os.path.join(temp_workspace_dir, 'taskgraph')
    for path in [results_dir, taskgraph_working_dir, temp_workspace_dir]:
        try:
            os.makedirs(results_dir)
        except OSError:
            pass
    task_graph = taskgraph.TaskGraph(taskgraph_working_dir, -1)

    aligned_pop_raster_list = align_raster_list(
        [x[0] for x in population_raster_path_id_list] + [dem_raster_path],
        temp_workspace_dir)

    for pop_index, (_, pop_id) in enumerate(population_raster_path_id_list):
        # mask to < 10m
        pop_height_masked_path = os.path.join(
            temp_workspace_dir, '%s_masked_by_10m.tif' % pop_id)
        raster_info = pygeoprocessing.get_raster_info(
            aligned_pop_raster_list[pop_index])

        pop_height_mask_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(aligned_pop_raster_list[pop_index], 1),
                 (aligned_pop_raster_list[-1], 1),
                 (10.0, 'raw'),  # mask to 10 meters
                 (raster_info['nodata'][0], 'raw')],  # the -1 index is the dem
                mask_by_height_op, pop_height_masked_path, gdal.GDT_Float32,
                raster_info['nodata'][0]),
            target_path_list=[pop_height_masked_path],
            task_name='pop height mask 10m %s' % pop_id)
        pop_height_mask_task.join()

        target_pixel_size = raster_info['pixel_size']
        # spread the < 10m population out 2km
        n_pixels_in_2km = int(2000.0 / (
            M_PER_DEGREE * abs(target_pixel_size[0])))
        kernel_radius_2km = [n_pixels_in_2km, n_pixels_in_2km]
        kernel_2km_filepath = os.path.join(
            temp_workspace_dir, '2km_kernel.tif')
        create_averaging_kernel_raster(
            kernel_radius_2km, kernel_2km_filepath, normalize=False)
        pop_sum_within_2km_path = os.path.join(
            temp_workspace_dir, '%s_pop_sum_within_2km.tif' % pop_id)
        pop_sum_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=(
                (pop_height_masked_path, 1), (kernel_2km_filepath, 1),
                pop_sum_within_2km_path),
            kwargs={'working_dir': temp_workspace_dir},
            target_path_list=[pop_sum_within_2km_path],
            task_name='pop sum w/in 2km %s' % pop_id)

        # spread the 2km pop out by the hab distance
        for habitat_id in habitat_fieldname_list:
            hab_raster_path, _, prot_distance = (
                habitat_vector_path_map[habitat_id])
            # make a kernel that goes out the distance of the protective
            # distance of habitat
            n_pixels_in_prot_dist = max(1, int(prot_distance / (
                M_PER_DEGREE * abs(target_pixel_size[0]))))
            kernel_radius = [n_pixels_in_prot_dist, n_pixels_in_prot_dist]
            kernel_filepath = os.path.join(
                temp_workspace_dir, '%s_kernel.tif' % habitat_id)
            create_averaging_kernel_raster(
                kernel_radius, kernel_filepath, normalize=False)
            population_hab_spread_raster_path = os.path.join(
                temp_workspace_dir, '%s_%s_spread.tif' % (habitat_id, pop_id))
            spread_to_hab_task = task_graph.add_task(
                func=clean_convolve_2d,
                args=(
                    (pop_sum_within_2km_path, 1), (kernel_filepath, 1),
                    population_hab_spread_raster_path),
                kwargs={'working_dir': temp_workspace_dir},
                target_path_list=[population_hab_spread_raster_path],
                dependent_task_list=[pop_sum_task],
                task_name='spread pop to hab %s %s ' % (pop_id, habitat_id))

            hab_raster_info = pygeoprocessing.get_raster_info(hab_raster_path)

            # warp pop result to overlay
            clipped_pop_hab_spread_raster_path = os.path.join(
                temp_workspace_dir, '%s_%s_spread_clipped.tif' % (
                    habitat_id, pop_id))
            task_graph.add_task(
                func=pygeoprocessing.warp_raster,
                args=(
                    population_hab_spread_raster_path,
                    hab_raster_info['pixel_size'],
                    clipped_pop_hab_spread_raster_path,
                    'near'),
                kwargs={
                    'target_bb': hab_raster_info['bounding_box'],
                    'working_dir': temp_workspace_dir},
                target_path_list=[clipped_pop_hab_spread_raster_path],
                dependent_task_list=[spread_to_hab_task],
                task_name='spread to hab %s %s' % (pop_id, habitat_id))

            # mask the convolution by the habitat mask
            task_graph.join()
            hab_spread_nodata = pygeoprocessing.get_raster_info(
                clipped_pop_hab_spread_raster_path)['nodata'][0]
            hab_nodata = pygeoprocessing.get_raster_info(
                hab_raster_path)['nodata'][0]
            habitat_value_raster_path = os.path.join(
                results_dir, '%s_%s_coverage.tif' % (habitat_id, pop_id))

            buffered_point_raster_mask_path = os.path.join(
                temp_workspace_dir, '%s_buffer_mask.tif' % habitat_id)
            buffered_raster_task = task_graph.add_task(
                func=make_buffered_point_raster_mask,
                args=(
                    shore_sample_point_vector_path,
                    clipped_pop_hab_spread_raster_path,
                    temp_workspace_dir, habitat_id,
                    prot_distance, buffered_point_raster_mask_path),
                target_path_list=[buffered_point_raster_mask_path],
                task_name='buffered point for %s' % habitat_id)

            hab_value_task = task_graph.add_task(
                func=pygeoprocessing.raster_calculator,
                args=(
                    [(clipped_pop_hab_spread_raster_path, 1),
                     (hab_raster_path, 1),
                     (buffered_point_raster_mask_path, 1),
                     (hab_spread_nodata, 'raw'),
                     (hab_nodata, 'raw')],
                    intersect_and_mask_raster_op, habitat_value_raster_path,
                    gdal.GDT_Float32, hab_spread_nodata),
                dependent_task_list=[buffered_raster_task],
                target_path_list=[habitat_value_raster_path],
                task_name='mask result %s %s' % (pop_id, habitat_id))

            task_graph.add_task(
                func=ecoshard.build_overviews,
                args=(habitat_value_raster_path,),
                target_path_list=[habitat_value_raster_path],
                dependent_task_list=[hab_value_task],
                task_name='build overviews for %s' % habitat_value_raster_path)

    task_graph.join()
    task_graph.close()
    del task_graph

def calculate_habitat_value(
        shore_sample_point_vector, template_raster_path,
        habitat_fieldname_list, habitat_vector_path_map, results_dir):
    """Calculate habitat value.

    Will create rasters in the `results_dir` directory named from the
    `habitat_fieldname_list` values containing relative importance of
    global habitat. The higher the value of a pixel the more important that
    pixel of habitat is for protection of the coastline.

    Parameters:
        shore_sample_point_vector (str): path to CV analysis vector containing
            at least the fields `Rt` and `Rt_nohab_[hab]` for all habitat
            types under consideration.
        template_raster_path (str): path to an existing raster whose size and
            shape will be used to be the base of the raster that's created
            for each habitat type.
        habitat_fieldname_list (list): list of habitat ids to analyise.
        habitat_vector_path_map (dict): maps fieldnames from
            `habitat_fieldname_list` to 3-tuples of
            (path to hab vector (str), risk val (float),
             protective distance (float)).
        results_dir (str): path to directory containing habitat back projection
            results

    Returns:
        None.

    """
    temp_workspace_dir = os.path.join(results_dir, 'hab_value_churn')
    for dir_path in [results_dir, temp_workspace_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    gpkg_driver = ogr.GetDriverByName('gpkg')
    shore_sample_point_vector = gdal.OpenEx(
        shore_sample_point_vector, gdal.OF_VECTOR)
    shore_sample_point_layer = shore_sample_point_vector.GetLayer()

    for habitat_id in habitat_fieldname_list:
        habitat_service_id = 'Rt_habservice_%s' % habitat_id
        hab_raster_path, _, protective_distance = (
            habitat_vector_path_map[habitat_id])

        buffer_habitat_path = os.path.join(
            temp_workspace_dir, '%s_buffer.gpkg' % habitat_id)
        buffer_habitat_vector = gpkg_driver.CreateDataSource(
            buffer_habitat_path)
        wgs84_srs = osr.SpatialReference()
        wgs84_srs.ImportFromEPSG(4326)
        buffer_habitat_layer = (
            buffer_habitat_vector.CreateLayer(
                habitat_service_id, wgs84_srs, ogr.wkbPolygon))
        buffer_habitat_layer.CreateField(ogr.FieldDefn(
            habitat_service_id, ogr.OFTReal))
        buffer_habitat_layer_defn = buffer_habitat_layer.GetLayerDefn()

        shore_sample_point_layer.ResetReading()
        buffer_habitat_layer.StartTransaction()
        for point_index, point_feature in enumerate(shore_sample_point_layer):
            if point_index % 1000 == 0:
                LOGGER.debug(
                    'point buffering is %.2f%% complete',
                    point_index / shore_sample_point_layer.GetFeatureCount() *
                    100.0)
            # for each point, convert to local UTM to buffer out a given
            # distance then back to wgs84
            point_geom = point_feature.GetGeometryRef()
            utm_srs = calculate_utm_srs(point_geom.GetX(), point_geom.GetY())
            wgs84_to_utm_transform = osr.CoordinateTransformation(
                wgs84_srs, utm_srs)
            utm_to_wgs84_transform = osr.CoordinateTransformation(
                utm_srs, wgs84_srs)
            point_geom.Transform(wgs84_to_utm_transform)
            buffer_poly_geom = point_geom.Buffer(protective_distance)
            buffer_poly_geom.Transform(utm_to_wgs84_transform)

            buffer_point_feature = ogr.Feature(buffer_habitat_layer_defn)
            buffer_point_feature.SetGeometry(buffer_poly_geom)

            # reefs are special
            if habitat_id in REEF_FIELDS:
                buffer_point_feature.SetField(
                    habitat_service_id,
                    point_feature.GetField('Rt_habservice_reefs_all'))
            else:
                buffer_point_feature.SetField(
                    habitat_service_id,
                    point_feature.GetField(habitat_service_id))
            buffer_habitat_layer.CreateFeature(buffer_point_feature)
            buffer_point_feature = None
            point_feature = None
            buffer_poly_geom = None
            point_geom = None

        # at this point every shore point has been buffered to the effective
        # habitat distance and the habitat service has been saved with it
        buffer_habitat_layer.CommitTransaction()
        buffer_habitat_layer = None
        buffer_habitat_vector = None
        value_coverage_raster_path = os.path.join(
            temp_workspace_dir, '%s_value_cover.tif' % habitat_id)
        pygeoprocessing.new_raster_from_base(
            template_raster_path, value_coverage_raster_path,
            gdal.GDT_Float32, [0],
            raster_driver_creation_tuple=(
                'GTIFF', (
                    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
                    'BLOCKXSIZE=256', 'BLOCKYSIZE=256', 'SPARSE_OK=TRUE')))
        pygeoprocessing.rasterize(
            buffer_habitat_path, value_coverage_raster_path,
            option_list=[
                'ATTRIBUTE=%s' % habitat_service_id,
                'MERGE_ALG=ADD'])

        habitat_value_raster_path = os.path.join(
            results_dir, '%s_value.tif' % habitat_id)

        value_coverage_nodata = pygeoprocessing.get_raster_info(
            value_coverage_raster_path)['nodata'][0]
        hab_nodata = pygeoprocessing.get_raster_info(
            hab_raster_path)['nodata'][0]

        aligned_value_hab_raster_path_list = align_raster_list(
            [value_coverage_raster_path, hab_raster_path],
            temp_workspace_dir)

        pygeoprocessing.raster_calculator(
            [(aligned_value_hab_raster_path_list[0], 1),
             (aligned_value_hab_raster_path_list[1], 1),
             (value_coverage_nodata, 'raw'), (hab_nodata, 'raw')],
            intersect_raster_op, habitat_value_raster_path, gdal.GDT_Float32,
            value_coverage_nodata)

        ecoshard.build_overviews(habitat_value_raster_path)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calcualte risk from reefs')
    parser.add_argument(
        'cv_risk_vector_pattern', nargs='+',
        help='Can be a pattern to a file.')
    args = parser.parse_args()

    for dir_path in [WORKSPACE_DIR, ECOSHARD_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, -1)
    tdd_downloader = taskgraph_downloader_pnn.TaskGraphDownloader(
        ECOSHARD_DIR, task_graph)

    tdd_downloader.download_ecoshard(
        GLOBAL_REEFS_RASTER_URL, 'reefs')
    tdd_downloader.download_ecoshard(
        LS_POPULATION_RASTER_URL, 'ls_population')
    tdd_downloader.download_ecoshard(
        POVERTY_POPULATION_RASTER_URL, 'poor_population')

    for cv_risk_vector_pattern in args.cv_risk_vector_pattern:
        for cv_path in glob.glob(cv_risk_vector_pattern):
            basename = os.path.basename(os.path.splitext(cv_path)[0])
            LOGGER.debug(basename)

            reefs_value_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_value_%s.tif" % basename)
            reefs_poor_pop_coverage_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_poor_pop_coverage_%s.tif" % basename)
            reefs_total_pop_coverage_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_total_pop_coverage_%s.tif" % basename)

    task_graph.join()
    task_graph.close()
