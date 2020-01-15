"""Global CV Analysis for CNC.

Design doc is here:

https://docs.google.com/document/d/18AcJM-rXeIYgEsmqlaUwdtm7gdWLiaD6kkRkpARILlw/edit#heading=h.bbujb61ete53

"""
import collections
import gzip
import logging
import math
import multiprocessing
import os
import shutil
import sys
import tempfile
import threading
import zipfile

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import numpy
import pygeoprocessing
import retrying
import rtree
import shapely.geometry
import shapely.strtree
import shapely.wkt
import taskgraph

WORKSPACE_DIR = 'global_cv_workspace'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
GRID_WORKSPACE_DIR = os.path.join(WORKSPACE_DIR, 'grid_workspaces')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
TARGET_NODATA = -1
TARGET_CV_VECTOR_PATH = os.path.join(
    WORKSPACE_DIR, 'global_cv_analysis_result.gpkg')

# [minx, miny, maxx, maxy].
GLOBAL_AOI_WGS84_BB = [-179, -65, 180, 77]
SHORE_POINT_SAMPLE_DISTANCE = 2000.0
RELIEF_SAMPLE_DISTANCE = 5000.0
N_FETCH_RAYS = 16  # this is hardcoded because of WWIII fields
MAX_FETCH_DISTANCE = 60000

ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')

GLOBAL_POLYGON_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_global_polygon_simplified_geometries_'
    'md5_653118dde775057e24de52542b01eaee.gpkg')

LULC_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_'
    'md5_1254d25f937e6d9bdee5779d377c5aa4.tif')
BUFFER_VECTOR_URL = (
    ECOSHARD_BUCKET_URL +
    'buffered_global_shore_5km_md5_a68e1049c1c03673add014cd29b7b368.gpkg')
SHORE_GRID_URL = (
    ECOSHARD_BUCKET_URL +
    'shore_grid_md5_1a09b69c0c548d3b25d7e14c3ddb60c9.gpkg')

GLOBAL_WWIII_GZ_URL = (
    ECOSHARD_BUCKET_URL +
    'wave_watch_iii_md5_c8bb1ce4739e0a27ee608303c217ab5b.gpkg.gz')
GLOBAL_DEM_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'global_dem_md5_22c5c09ac4c4c722c844ab331b34996c.tif')
LS_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')
SLR_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'MSL_Map_MERGED_Global_AVISO_NoGIA_Adjust_'
    'md5_3072845759841d0b2523d00fe9518fee.tif')
GLOBAL_GEOMORPHOLOGY_VECTOR_URL = (
    ECOSHARD_BUCKET_URL +
    'geomorphology_md5_ace4406c75459c56f158c9f5a4636897.gpkg')
GLOBAL_REEFS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_reef_md5_5a90d55a505813b5aa9662faee351bf8.tif')
GLOBAL_MANGROVES_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_mangrove_md5_0ec85cb51dab3c9ec3215783268111cc.tif')
GLOBAL_SEAGRASS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_seagrass_md5_a9cc6d922d2e74a14f74b4107c94a0d6.tif')
GLOBAL_SALTMARSH_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_saltmarsh_md5_203d8600fd4b6df91f53f66f2a011bcd.tif')

GLOBAL_DATA_URL_MAP = {
    'geomorphology': GLOBAL_GEOMORPHOLOGY_VECTOR_URL,
    'mangroves': GLOBAL_MANGROVES_RASTER_URL,
    'reefs': GLOBAL_REEFS_RASTER_URL,
    'seagrass': GLOBAL_SEAGRASS_RASTER_URL,
    'saltmarsh': GLOBAL_SALTMARSH_RASTER_URL,
    'dem': GLOBAL_DEM_RASTER_URL,
    'slr': SLR_RASTER_URL,
    'landmass': GLOBAL_POLYGON_URL,
    'shore_grid': SHORE_GRID_URL,
    'lulc': LULC_RASTER_URL
    }

SEDTYPE_TO_RISK = {
    0: 5,  # unknown
    1: 5,  # sandy
    2: 1,  # unerodable
    3: 5,  # muddy
    4: 1,  # coral/mangrove
}

# This dictionary maps landcode id to (risk, dist) tuples
LULC_CODE_TO_HAB_MAP = {
    0: (0, None),
    10: (0, None),
    11: (0, None),
    12: (0, None),
    20: (0, None),
    30: (0, None),
    40: (0, None),
    50: (1, 2000),
    60: (1, 2000),
    61: (1, 2000),
    62: (1, 2000),
    70: (1, 2000),
    71: (1, 2000),
    72: (1, 2000),
    80: (1, 2000),
    81: (1, 2000),
    82: (1, 2000),
    90: (1, 2000),
    100: (1, 2000),
    110: (2, 2000),
    120: (2, 2000),
    121: (2, 2000),
    122: (2, 2000),
    130: (2, 2000),
    140: (2, 2000),
    150: (4, 500),
    151: (4, 500),
    152: (4, 500),
    153: (4, 500),
    160: (2, 1000),
    170: (2, 1000),
    180: (2, 1000),
    190: (0, None),
    200: (0, None),
    201: (0, None),
    202: (0, None),
    210: (0, None),
    220: (0, None),
    }

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

STOP_SENTINEL = 'STOP'

HABITAT_VECTOR_PATH_MAP = {
    'reefs': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['reefs'])),
        1, 2000.0),
    'mangroves': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['mangroves'])),
        1, 1000.0),
    'saltmarsh': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['saltmarsh'])),
        2, 1000.0),
    'seagrass': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['seagrass'])),
        4, 500.0)}


def download_and_unzip(url, target_dir, target_token_path):
    """Download `url` to `target_dir` and touch `target_token_path`."""
    zipfile_path = os.path.join(target_dir, os.path.basename(url))
    LOGGER.debug('url %s, zipfile_path: %s', url, zipfile_path)
    ecoshard.download_url(url, zipfile_path)

    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(target_token_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def download_and_ungzip(url, target_path, buffer_size=2**20):
    """Download `url` to `target_dir` and touch `target_token_path`."""
    gzipfile_path = os.path.join(
        os.path.dirname(target_path), os.path.basename(url))
    ecoshard.download_url(url, gzipfile_path)

    with gzip.open(gzipfile_path, 'rb') as gzip_file:
        with open(target_path, 'wb') as target_file:
            while True:
                content = gzip_file.read(buffer_size)
                if content:
                    target_file.write(content)
                else:
                    break


def build_rtree(vector_path):
    """Build an rtree that can be queried for nearest neighbors.

    Parameters:
        vector_path (str): path to vector of geometry to build into
            r tree.

    Returns:
        rtree.Index object that will return shapely geometry objects with
            a field_val_map field that contains the 'fieldname'->value pairs
            from the original vector. The main object will also have a
            `field_name_type_list` field which contains original
            fieldname/field type pairs

    """
    geometry_prep_list = []
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_name_type_list = []
    for index in range(layer_defn.GetFieldCount()):
        field_name = layer_defn.GetFieldDefn(index).GetName()
        field_type = layer_defn.GetFieldDefn(index).GetType()
        field_name_type_list.append((field_name, field_type))

    LOGGER.debug('loop through features for rtree')
    for index, feature in enumerate(layer):
        feature_geom = feature.GetGeometryRef().Clone()
        feature_geom_shapely = shapely.wkb.loads(feature_geom.ExportToWkb())
        field_val_map = {}
        for field_name, _ in field_name_type_list:
            field_val_map[field_name] = (
                feature.GetField(field_name))
        geometry_prep_list.append(
            (index, feature_geom_shapely.bounds, field_val_map))
    LOGGER.debug('constructing the tree')
    r_tree = rtree.index.Index(geometry_prep_list)
    LOGGER.debug('all done')
    r_tree.field_name_type_list = field_name_type_list
    return r_tree


def build_strtree(vector_path):
    """Build an rtree that generates geom and preped geometry.

    Parameters:
        vector_path (str): path to vector of geometry to build into
            r tree.

    Returns:
        strtree.STRtree object that will return shapely geometry objects
            with a .prep field that is prepared geomtry for fast testing,
            a .geom field that is the base gdal geometry, and a field_val_map
            field that contains the 'fieldname'->value pairs from the original
            vector. The main object will also have a `field_name_type_list`
            field which contains original fieldname/field type pairs

    """
    geometry_prep_list = []
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_name_type_list = []
    for index in range(layer_defn.GetFieldCount()):
        field_name = layer_defn.GetFieldDefn(index).GetName()
        field_type = layer_defn.GetFieldDefn(index).GetType()
        field_name_type_list.append((field_name, field_type))

    LOGGER.debug('loop through features for rtree')
    for index, feature in enumerate(layer):
        feature_geom = feature.GetGeometryRef().Clone()
        feature_geom_shapely = shapely.wkb.loads(feature_geom.ExportToWkb())
        feature_geom_shapely.prep = shapely.prepared.prep(feature_geom_shapely)
        feature_geom_shapely.geom = feature_geom
        feature_geom_shapely.id = index
        feature_geom_shapely.field_val_map = {}
        for field_name, _ in field_name_type_list:
            feature_geom_shapely.field_val_map[field_name] = (
                feature.GetField(field_name))
        geometry_prep_list.append(feature_geom_shapely)
    LOGGER.debug('constructing the tree')
    r_tree = shapely.strtree.STRtree(geometry_prep_list)
    LOGGER.debug('all done')
    r_tree.field_name_type_list = field_name_type_list
    return r_tree


def cv_grid_worker(
        bb_work_queue,
        cv_point_complete_queue,
        global_landmass_vector_path,
        geomorphology_vector_path,
        slr_raster_path,
        global_dem_raster_path,
        wwiii_vector_path,
        habitat_raster_path_map,
        ):
    """Worker process to calculate CV for a grid.

    Parameters:
        bb_work_queue (multiprocessing.Queue): contains
            [minx, miny, maxx, maxy] bounding box values to be processed or
            `STOP_SENTINEL` values to indicate the worker should be terminated.
        cv_point_complete_queue (multiprocessing.Queue): this queue is used to
            pass completed CV point vectors for further processing. It will be
            terminated by `STOP_SENTINEL`.
        global_landmass_vector_path (str): path to global landmass vector. Used
            for intersecting rays to determine shoreline exposure.
        geomorphology_vector_path (str): path to line geometry geomorphology
            layer that has a field called 'SEDTYPE'
        global_dem_raster_path (str): path to a global dem/bathymetry raster.
        slr_raster_path (str): path to a sea level rise raster.
        wwiii_vector_path (str): path to wave watch III dataset that has
            fields TODO FILL IN WHAT THEY ARE
        habitat_raster_path_map (dict): mapt a habitat id to a
            (raster path, risk, dist(m)) tuple. These are the raster versions
            of habitats to use in Rhab.

    Returns:
        None.

    """
    LOGGER.info('build geomorphology rtree')
    geomorphology_strtree = build_strtree(geomorphology_vector_path)
    geomorphology_proj_wkt = pygeoprocessing.get_vector_info(
        geomorphology_vector_path)['projection']
    gegeomorphology_proj = osr.SpatialReference()
    gegeomorphology_proj.ImportFromWkt(geomorphology_proj_wkt)

    LOGGER.info('build landmass rtree')
    landmass_strtree = build_strtree(global_landmass_vector_path)

    LOGGER.info('build wwiii rtree')
    wwiii_rtree = build_rtree(wwiii_vector_path)

    target_pixel_size = [
        SHORE_POINT_SAMPLE_DISTANCE / 4,
        -SHORE_POINT_SAMPLE_DISTANCE / 4]

    while True:
        payload = bb_work_queue.get()
        if payload == STOP_SENTINEL:
            LOGGER.debug('stopping')
            # put it back so others can stop
            bb_work_queue.put(STOP_SENTINEL)
            break
        else:
            LOGGER.debug('running')
        # otherwise payload is the bounding box
        index, (lng_min, lat_min, lng_max, lat_max) = payload
        bounding_box_list = [lng_min, lat_min, lng_max, lat_max]
        # create workspace
        workspace_dir = os.path.join(
            GRID_WORKSPACE_DIR, '%d_%s_%s_%s_%s' % (
                index, lng_min, lat_min, lng_max, lat_max))

        try:
            os.makedirs(workspace_dir)
        except OSError:
            pass

        # task_graph = taskgraph.TaskGraph(workspace_dir, -1)

        utm_srs = calculate_utm_srs((lng_min+lng_max)/2, (lat_min+lat_max)/2)
        wgs84_srs = osr.SpatialReference()
        wgs84_srs.ImportFromEPSG(4326)

        try:
            local_geomorphology_vector_path = os.path.join(
                workspace_dir, 'geomorphology.gpkg')
            clip_geometry(
                bounding_box_list, wgs84_srs, utm_srs,
                ogr.wkbMultiLineString, geomorphology_strtree,
                local_geomorphology_vector_path)

            local_landmass_vector_path = os.path.join(
                workspace_dir, 'landmass.gpkg')
            clip_geometry(
                bounding_box_list, wgs84_srs, utm_srs,
                ogr.wkbPolygon, landmass_strtree,
                local_landmass_vector_path)

            landmass_boundary_vector_path = os.path.join(
                workspace_dir, 'landmass_boundary.gpkg')
            vector_to_lines(
                local_landmass_vector_path, landmass_boundary_vector_path)

            local_dem_path = os.path.join(
                workspace_dir, 'dem.tif')
            clip_and_reproject_raster(
                global_dem_raster_path, local_dem_path, utm_srs.ExportToWkt(),
                bounding_box_list, RELIEF_SAMPLE_DISTANCE, 'bilinear', True,
                target_pixel_size)

            local_slr_path = os.path.join(
                workspace_dir, 'slr.tif')
            clip_and_reproject_raster(
                global_dem_raster_path, local_slr_path, utm_srs.ExportToWkt(),
                bounding_box_list, 0, 'bilinear', True,
                target_pixel_size)

            shore_point_vector_path = os.path.join(
                workspace_dir, 'shore_points.gpkg')
            sample_line_to_points(
                local_geomorphology_vector_path, shore_point_vector_path,
                SHORE_POINT_SAMPLE_DISTANCE)

            # Rrelief
            LOGGER.info('calculate relief on %s', workspace_dir)
            calculate_relief(
                shore_point_vector_path, local_dem_path, 'relief')
            LOGGER.info('calculate rhab on %s', workspace_dir)
            # Rhab
            calculate_rhab(
                shore_point_vector_path, habitat_raster_path_map, 'Rhab',
                target_pixel_size)

            # Rslr
            calculate_slr(shore_point_vector_path, local_slr_path, 'slr')

            # wind and wave power
            calculate_wind_and_wave(
                shore_point_vector_path, landmass_boundary_vector_path,
                local_dem_path, wwiii_rtree, 'rei', 'ew')

            # Rsurge
            calculate_surge(shore_point_vector_path, local_dem_path, 'surge')

            # Rgeomorphology
            calculate_geomorphology(
                shore_point_vector_path, local_geomorphology_vector_path,
                'Rgeomorphology')

            LOGGER.info('completed %s', shore_point_vector_path)
            cv_point_complete_queue.put(shore_point_vector_path)

        except ValueError as e:
            if 'no data intersects this box' in str(e):
                LOGGER.exception('error on %s', payload)
                LOGGER.warning('missing data, removing workspace')
                retrying_rmtree(workspace_dir)
            else:
                raise


def make_shore_kernel(kernel_path):
    """Make a 3x3 raster with a 9 in the middle and 1s on the outside."""
    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_path.encode('utf-8'), 3, 3, 1,
        gdal.GDT_Byte)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([0, 1, 0, 0, 0, -1])
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(127)
    kernel_band.WriteArray(numpy.array([[1, 1, 1], [1, 9, 1], [1, 1, 1]]))


def calculate_geomorphology(
        shore_point_vector_path, geomorphology_vector_path,
        geomorphology_fieldname):
    """Sample the geomorphology vector path for the closest line to each point.

    Parameters:
        shore_point_vector_path (str):  path to a point shapefile to
            for relief point analysis.
        geomorphology_vector_path (str): path to a vector of lines that
            contains the integer field 'SEDTYPE'.
        geomorphology_fieldname (str): fieldname to add to
            `shore_point_vector_path`.

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()
    shore_point_layer.CreateField(
        ogr.FieldDefn(geomorphology_fieldname, ogr.OFTReal))
    geomorphology_strtree = build_strtree(geomorphology_vector_path)

    shore_point_layer.StartTransaction()
    for shore_point_feature in shore_point_layer:
        shore_point_geom = shapely.wkb.loads(
            shore_point_feature.GetGeometryRef().ExportToWkb())
        min_dist = MAX_FETCH_DISTANCE
        geo_risk = 5
        for line in geomorphology_strtree.query(shore_point_geom.buffer(500)):
            cur_dist = line.distance(shore_point_geom)
            if cur_dist < min_dist:
                min_dist = cur_dist
                geo_risk = line.field_val_map['Rgeo']
        shore_point_feature.SetField(geomorphology_fieldname, geo_risk)
        shore_point_layer.SetFeature(shore_point_feature)
    shore_point_layer.CommitTransaction()
    shore_point_layer = None
    shore_point_vector = None


def calculate_surge(
        shore_point_vector_path, bathymetry_raster_path, surge_fieldname):
    """Calculate surge potential as distance to continental shelf (-150m).

    Parameters:
        base_shore_point_vector_path (string):  path to a point shapefile to
            for relief point analysis.
        global_dem_path (string): path to a DEM raster projected in wgs84.
        surge_fieldname (str): fieldname to add to `shore_point_vector_path`
        workspace_dir (string): path to a directory to make local calculations
            in
        target_surge_point_vector_path (string): path to output vector.
            after completion will a value for closest distance to continental
            shelf called 'surge'.

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()
    shore_point_layer.CreateField(ogr.FieldDefn(surge_fieldname, ogr.OFTReal))

    shelf_nodata = 2

    bathymetry_nodata = pygeoprocessing.get_raster_info(
        bathymetry_raster_path)['nodata'][0]

    def mask_shelf(depth_array):
        valid_mask = ~numpy.isclose(depth_array, bathymetry_nodata)
        result_array = numpy.empty(
            depth_array.shape, dtype=numpy.int16)
        result_array[:] = shelf_nodata
        result_array[valid_mask] = 0
        result_array[depth_array < -150] = 1
        return result_array

    workspace_dir = os.path.dirname(shore_point_vector_path)

    shelf_mask_path = os.path.join(workspace_dir, 'shelf_mask.tif')
    pygeoprocessing.raster_calculator(
        [(bathymetry_raster_path, 1)], mask_shelf,
        shelf_mask_path, gdal.GDT_Byte, shelf_nodata)

    # convolve to find edges
    # grid shoreline from raster
    shelf_kernel_path = os.path.join(workspace_dir, 'shelf_kernel.tif')
    shelf_convoultion_raster_path = os.path.join(
        workspace_dir, 'shelf_convolution.tif')
    make_shore_kernel(shelf_kernel_path)
    pygeoprocessing.convolve_2d(
        (shelf_mask_path, 1), (shelf_kernel_path, 1),
        shelf_convoultion_raster_path, target_datatype=gdal.GDT_Byte,
        target_nodata=255)

    nodata = pygeoprocessing.get_raster_info(
        shelf_convoultion_raster_path)['nodata'][0]

    def _shelf_mask_op(shelf_convolution):
        """Mask values on land that border the continental shelf."""
        result = numpy.empty(shelf_convolution.shape, dtype=numpy.uint8)
        result[:] = nodata
        valid_mask = shelf_convolution != nodata
        # If a pixel is on land, it gets at least a 9, but if it's all on
        # land it gets an 17 (8 neighboring pixels), so we search between 9
        # and 17 to determine a shore pixel
        result[valid_mask] = numpy.where(
            (shelf_convolution[valid_mask] >= 9) &
            (shelf_convolution[valid_mask] < 17), 1, nodata)
        return result

    shelf_edge_raster_path = os.path.join(workspace_dir, 'shelf_edge.tif')
    pygeoprocessing.raster_calculator(
        [(shelf_convoultion_raster_path, 1)], _shelf_mask_op,
        shelf_edge_raster_path, gdal.GDT_Byte, nodata)

    shore_geotransform = pygeoprocessing.get_raster_info(
        shelf_edge_raster_path)['geotransform']

    shelf_rtree = rtree.index.Index()

    for offset_info, data_block in pygeoprocessing.iterblocks(
            (shelf_edge_raster_path, 1)):
        row_indexes, col_indexes = numpy.mgrid[
            offset_info['yoff']:offset_info['yoff']+offset_info['win_ysize'],
            offset_info['xoff']:offset_info['xoff']+offset_info['win_xsize']]
        valid_mask = data_block == 1
        x_coordinates = (
            shore_geotransform[0] +
            shore_geotransform[1] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[2] * (row_indexes[valid_mask] + 0.5))
        y_coordinates = (
            shore_geotransform[3] +
            shore_geotransform[4] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[5] * (row_indexes[valid_mask] + 0.5))

        for x_coord, y_coord in zip(x_coordinates, y_coordinates):
            shelf_rtree.insert(
                0, [x_coord, y_coord, x_coord, y_coord],
                obj=shapely.geometry.Point(x_coord, y_coord))

    shore_point_layer.StartTransaction()
    for point_feature in shore_point_layer:
        point_geometry = point_feature.GetGeometryRef()
        point_shapely = shapely.wkb.loads(point_geometry.ExportToWkb())
        nearest_point = list(shelf_rtree.nearest(
                (point_geometry.GetX(),
                 point_geometry.GetY(),
                 point_geometry.GetX(),
                 point_geometry.GetY()),
                objects='raw', num_results=1))
        if len(nearest_point) > 0:
            distance = nearest_point[0].distance(point_shapely)
            point_feature.SetField(surge_fieldname, float(distance))
        else:
            # so far away it's essentially not an issue
            point_feature.SetField(surge_fieldname, MAX_FETCH_DISTANCE)
        shore_point_layer.SetFeature(point_feature)

    shore_point_layer.CommitTransaction()
    shore_point_layer.SyncToDisk()
    shore_point_layer = None
    shore_point_vector = None


def calculate_wind_and_wave(
        shore_point_vector_path, landmass_vector_path, bathymetry_raster_path,
        wwiii_rtree, wind_fieldname, wave_fieldname):
    """Calculate wind exposure for given points.

    Parameters:
        shore_point_vector_path (str): path to a point vector, this value will
            be modified to hold the total wind exposure at this point.
        landmass_vector_path (str): path to a vector indicating landmass that
            will block wind exposure.
        bathymetry_raster_path (str): path to a raster indicating bathymetry
            values. (negative is deeper).
        wwiii_rtree (str): path to an r_tree that can find the nearest point
            in lat/lng whose object has values 'REI_PCT', 'REI_V',
            'WavP_[DIR]', 'WavPPCT', 'V10PCT_[DIR]'.
        wind_fieldname (str): fieldname to add to `shore_point_vector_path` for
            wind power.
        wave_fieldname (str): fieldname to add to `shore_point_vector_path` for
            wave power.

    Returns:
        None

    """
    gpkg_driver = ogr.GetDriverByName('gpkg')
    temp_workspace_dir = tempfile.mkdtemp(
        dir=os.path.dirname(shore_point_vector_path),
        prefix='calculate_rwind_')
    temp_fetch_rays_vector = gpkg_driver.CreateDataSource(
        os.path.join(temp_workspace_dir, 'fetch_rays.gpkg'))
    layer_name = 'fetch_rays'
    shore_point_projection_wkt = pygeoprocessing.get_vector_info(
        shore_point_vector_path)['projection']
    shore_point_srs = osr.SpatialReference()
    shore_point_srs.ImportFromWkt(shore_point_projection_wkt)
    temp_fetch_rays_layer = (
        temp_fetch_rays_vector.CreateLayer(
            str(layer_name), shore_point_srs, ogr.wkbLineString))
    temp_fetch_rays_layer.CreateField(ogr.FieldDefn(
        'fetch_dist', ogr.OFTReal))
    temp_fetch_rays_layer.CreateField(ogr.FieldDefn(
        'direction', ogr.OFTReal))
    temp_fetch_rays_defn = temp_fetch_rays_layer.GetLayerDefn()

    # These WWIII fields are the only ones needed for wind & wave equations
    # Copy them to a new vector which also gets more fields added with
    # computed values.
    target_shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    target_shore_point_layer = target_shore_point_vector.GetLayer()
    target_shore_point_layer.CreateField(
        ogr.FieldDefn(wind_fieldname, ogr.OFTReal))
    target_shore_point_layer.CreateField(
        ogr.FieldDefn(wave_fieldname, ogr.OFTReal))
    for ray_index in range(N_FETCH_RAYS):
        compass_degree = int(ray_index * 360 / N_FETCH_RAYS)
        target_shore_point_layer.CreateField(
            ogr.FieldDefn('fdist_%d' % compass_degree, ogr.OFTReal))
        target_shore_point_layer.CreateField(
            ogr.FieldDefn('fdepth_%d' % compass_degree, ogr.OFTReal))

    # Iterate over every shore point
    LOGGER.info("Casting rays and extracting bathymetry values")
    bathy_raster = gdal.OpenEx(
        bathymetry_raster_path, gdal.OF_RASTER | gdal.GA_ReadOnly)
    bathy_band = bathy_raster.GetRasterBand(1)
    bathy_raster_info = pygeoprocessing.get_raster_info(bathymetry_raster_path)
    bathy_gt = bathy_raster_info['geotransform']
    bathy_inv_gt = gdal.InvGeoTransform(bathy_gt)

    landmass_vector = gdal.OpenEx(landmass_vector_path, gdal.OF_VECTOR)
    landmass_layer = landmass_vector.GetLayer()
    landmass_geom_list = [
        shapely.wkb.loads(f.GetGeometryRef().ExportToWkb())
        for f in landmass_layer]
    landmass_union_geom = shapely.ops.cascaded_union(landmass_geom_list)
    landmass_layer = None
    landmass_vector = None
    landmass_union_geom_prep = shapely.prepared.prep(landmass_union_geom)

    landmass_strtree = build_strtree(landmass_vector_path)

    target_shore_point_layer.StartTransaction()
    temp_fetch_rays_layer.StartTransaction()

    # make a transfomer for local points to lat/lng for wwiii_rtree
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    base_to_target_transform = osr.CoordinateTransformation(
        shore_point_srs, wgs84_srs)

    for shore_point_feature in target_shore_point_layer:
        shore_point_geom = shore_point_feature.GetGeometryRef().Clone()
        _ = shore_point_geom.Transform(base_to_target_transform)
        wwiii_point = next(wwiii_rtree.nearest(
            (shore_point_geom.GetX(), shore_point_geom.GetY()), 1,
            objects='raw'))
        rei_value = 0.0
        height_list = []
        period_list = []
        e_local = 0.0
        e_ocean = 0.0
        # Iterate over every ray direction
        for sample_index in range(N_FETCH_RAYS):
            compass_degree = int(sample_index * 360 / N_FETCH_RAYS)
            compass_theta = float(sample_index) / N_FETCH_RAYS * 360

            # wwiii_point shoudl be closest point to shore point

            rei_pct = wwiii_point[
                'REI_PCT%d' % int(compass_theta)]
            rei_v = wwiii_point[
                'REI_V%d' % int(compass_theta)]
            cartesian_theta = -(compass_theta - 90)

            # Determine the direction the ray will point
            delta_x = math.cos(cartesian_theta * math.pi / 180)
            delta_y = math.sin(cartesian_theta * math.pi / 180)

            # Start a ray offset from the shore point
            # so that rays start outside of the landmass.
            # Shore points are interpolated onto the coastline,
            # but floating point error results in points being just
            # barely inside/outside the landmass.
            offset = 1  # 1 meter should be plenty
            shore_point_geometry = shore_point_feature.GetGeometryRef()
            point_a_x = (
                shore_point_geometry.GetX() + delta_x * offset)
            point_a_y = (
                shore_point_geometry.GetY() + delta_y * offset)
            point_b_x = point_a_x + delta_x * (
                MAX_FETCH_DISTANCE)
            point_b_y = point_a_y + delta_y * (
                MAX_FETCH_DISTANCE)
            shore_point_geometry = None

            # build ray geometry so we can intersect it later
            ray_geometry = ogr.Geometry(ogr.wkbLineString)
            ray_geometry.AddPoint(point_a_x, point_a_y)
            ray_geometry.AddPoint(point_b_x, point_b_y)

            # keep a shapely version of the ray so we can do fast intersection
            # with it and the entire landmass
            ray_point_origin_shapely = shapely.geometry.Point(
                point_a_x, point_a_y)

            if not landmass_union_geom_prep.intersects(
                    ray_point_origin_shapely):
                # the origin is in ocean, so we'll get a ray length > 0.0

                # This algorithm searches for intersections, if one is found
                # the ray updates and a smaller intersection set is determined
                # by experimentation I've found this is significant, but not
                # an order of magnitude, faster than looping through all
                # original possible intersections.  Since this algorithm
                # will be run for a long time, it's worth the additional
                # complexity
                tested_indexes = set()
                while True:
                    intersection = False
                    ray_envelope = ray_geometry.GetEnvelope()
                    for landmass_line in landmass_strtree.query(
                             shapely.geometry.box(
                                *[ray_envelope[i] for i in [0, 2, 1, 3]])):
                        if landmass_line.id in tested_indexes:
                            continue
                        tested_indexes.add(landmass_line.id)
                        if ray_geometry.Intersects(landmass_line.geom):
                            intersection_point = ray_geometry.Intersection(
                                landmass_line.geom)
                            # offset the dist with smallest_feature_size
                            # update the endpoint of the ray
                            ray_geometry = ogr.Geometry(ogr.wkbLineString)
                            ray_geometry.AddPoint(point_a_x, point_a_y)
                            ray_geometry.AddPoint(
                                intersection_point.GetX(),
                                intersection_point.GetY())
                            intersection = True
                            break
                    if not intersection:
                        break

                ray_step_loc = 0.0
                bathy_values = []
                # walk along ray
                ray_shapely = shapely.wkb.loads(ray_geometry.ExportToWkb())
                while ray_step_loc < ray_shapely.length:
                    sample_point = ray_shapely.interpolate(ray_step_loc)
                    ray_step_loc += SHORE_POINT_SAMPLE_DISTANCE/4
                    pixel_x, pixel_y = [int(x) for x in gdal.ApplyGeoTransform(
                        bathy_inv_gt,
                        sample_point.coords[0][0], sample_point.coords[0][1])]
                    if (pixel_x < 0 or pixel_y < 0 or
                            pixel_x >= bathy_band.XSize or
                            pixel_y >= bathy_band.YSize):
                        continue
                    bathy_values.append(
                        bathy_band.ReadAsArray(
                            pixel_x, pixel_y, 1, 1)[0][0])

                if bathy_values:
                    avg_bathy_value = numpy.mean(bathy_values)
                else:
                    avg_bathy_value = 0.0
                # when we get here, we have the final ray geometry
                ray_feature = ogr.Feature(temp_fetch_rays_defn)
                ray_feature.SetField('fetch_dist', ray_shapely.length)
                ray_feature.SetField('direction', compass_degree)
                ray_feature.SetGeometry(ray_geometry)
                temp_fetch_rays_layer.CreateFeature(ray_feature)
                rei_value += ray_shapely.length * rei_pct * rei_v
                ray_length = ray_geometry.Length()
                ray_feature = None
                ray_geometry = None

                shore_point_feature.SetField(
                    'fdist_%d' % compass_degree, ray_length)
                shore_point_feature.SetField(
                    'fdepth_%d' % compass_degree, float(avg_bathy_value))

                velocity = wwiii_point['V10PCT_%d' % compass_degree]
                occurrence = wwiii_point['REI_PCT%d' % compass_degree]

                height = compute_wave_height(
                    velocity, ray_shapely.length, avg_bathy_value)
                height_list.append(height)
                period = compute_wave_period(
                    velocity, ray_shapely.length, avg_bathy_value)
                period_list.append(period)
                power = 0.5 * float(height)**2 * float(period)  # UG Eq. 8
                e_local += power * occurrence  # UG Eq. 9

                if intersection:
                    e_ocean += (
                        wwiii_point['WavP_%d' % compass_degree] *
                        wwiii_point['WavPPCT%d' % compass_degree])

                ray_feature = None
                ray_geometry = None
                rei_value += ray_length * rei_pct * rei_v
            shore_point_feature.SetField(wind_fieldname, rei_value)
            shore_point_feature.SetField(wave_fieldname, max(e_ocean, e_local))
            target_shore_point_layer.SetFeature(shore_point_feature)

    target_shore_point_layer.CommitTransaction()
    target_shore_point_layer.SyncToDisk()
    target_shore_point_layer = None
    target_shore_point_vector = None
    temp_fetch_rays_layer.CommitTransaction()
    temp_fetch_rays_layer.SyncToDisk()
    temp_fetch_rays_layer = None
    temp_fetch_rays_vector = None
    bathy_raster = None
    bathy_band = None


def calculate_slr(shore_point_vector_path, slr_raster_path, target_fieldname):
    """Sample sea level rise raster and store values in shore points.

    Parameters:
        shore_point_vector_path (str): path to a vector of points in a local
            projected coordinate system. This vector will be modified by this
            function to include a new field called `target_fieldname`
            containing the weighted Rhab risk for the given point.
        slr_raster_path (str): path to a sea level rise raster indicating
            sea level rise amout in m.
        target_fieldname (str): fieldname to add to `shore_point_vector_path`
            that will contain the value of sea level rise for that point.

    Returns:
        None.

    """
    try:
        shore_point_vector = gdal.OpenEx(
            shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
        shore_point_layer = shore_point_vector.GetLayer()
        slr_field = ogr.FieldDefn(target_fieldname, ogr.OFTReal)
        slr_field.SetPrecision(5)
        slr_field.SetWidth(24)
        shore_point_layer.CreateField(slr_field)
        slr_info = pygeoprocessing.get_raster_info(slr_raster_path)
        inv_gt = gdal.InvGeoTransform(slr_info['geotransform'])

        slr_raster = gdal.OpenEx(slr_raster_path, gdal.OF_RASTER)
        slr_band = slr_raster.GetRasterBand(1)

        shore_point_layer.ResetReading()
        for point_feature in shore_point_layer:
            point_geometry = point_feature.GetGeometryRef()
            point_x, point_y = point_geometry.GetX(), point_geometry.GetY()
            point_geometry = None

            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(inv_gt, point_x, point_y)]

            try:
                pixel_value = slr_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception(
                    'slr_band size %d %d', slr_band.XSize,
                    slr_band.YSize)
                raise
            point_feature.SetField(target_fieldname, float(pixel_value))
            shore_point_layer.SetFeature(point_feature)

        shore_point_layer.SyncToDisk()
        shore_point_layer = None
        shore_point_vector = None
        slr_raster = None
        slr_band = None

    except Exception:
        LOGGER.exception('error in slr calc')
        raise


def calculate_rhab(
        shore_point_vector_path, habitat_raster_path_map, target_fieldname,
        target_pixel_size):
    """Add Rhab risk to the shore point vector path.

    Parameters:
        shore_point_vector_path (str): path to a vector of points in a local
            projected coordinate system. This vector will be modified by this
            function to include a new field called `target_fieldname`
            containing the weighted Rhab risk for the given point.
        habitat_raster_path_map (dict): a dictionary mapping "hab id"s to
            (path to raster, risk, effective distance) tuples.
        target_fieldname (str): fieldname to add to `shore_point_vector_path`
            that will contain the value of Rhab calculated for that point.
        target_pixel_size (list): x/y size of clipped habitat in projected
            units

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()
    relief_field = ogr.FieldDefn(target_fieldname, ogr.OFTReal)
    relief_field.SetPrecision(5)
    relief_field.SetWidth(24)
    shore_point_layer.CreateField(relief_field)

    shore_point_info = pygeoprocessing.get_vector_info(shore_point_vector_path)

    shore_point_feature_risk_map = collections.defaultdict(list)

    tmp_working_dir = tempfile.mkdtemp(
        prefix='calculate_rhab_',
        dir=os.path.dirname(shore_point_vector_path))
    LOGGER.debug(tmp_working_dir)
    for hab_id, (hab_raster_path, risk_val, eff_dist) in (
                habitat_raster_path_map.items()):
        local_hab_raster_path = os.path.join(
            tmp_working_dir, '%s.tif' % str(hab_id))
        LOGGER.debug(
            'clip %s to %s', hab_raster_path, shore_point_info['bounding_box'])
        clip_and_reproject_raster(
            hab_raster_path, local_hab_raster_path,
            shore_point_info['projection'],
            shore_point_info['bounding_box'], eff_dist, 'near', False,
            target_pixel_size)

        # make a convolution kernel as wide as the distance but adapted to
        # the non-square size of the rasters
        kernel_filepath = '%s_kernel%s' % os.path.splitext(
            local_hab_raster_path)
        kernel_radius = [abs(eff_dist / x) for x in target_pixel_size]
        create_averaging_kernel_raster(kernel_radius, kernel_filepath)
        hab_effective_area_raster_path = (
            '%s_effective_hab%s' % os.path.splitext(local_hab_raster_path))
        pygeoprocessing.convolve_2d(
            (local_hab_raster_path, 1), (kernel_filepath, 1),
            hab_effective_area_raster_path, mask_nodata=False)
        gt = pygeoprocessing.get_raster_info(
            hab_effective_area_raster_path)['geotransform']
        inv_gt = gdal.InvGeoTransform(gt)

        hab_effective_raster = gdal.OpenEx(
            hab_effective_area_raster_path, gdal.OF_RASTER)
        hab_effective_band = hab_effective_raster.GetRasterBand(1)
        shore_point_layer.ResetReading()
        for shore_feature in shore_point_layer:
            shore_geom = shore_feature.GetGeometryRef()
            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(
                    inv_gt, shore_geom.GetX(), shore_geom.GetY())]
            try:
                pixel_val = hab_effective_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception('error on pixel fetch for hab')
            if numpy.isclose(pixel_val, 0.0):
                pixel_val = 0
            # use max risk if no coverage
            shore_point_feature_risk_map[shore_feature.GetFID()].append(
                risk_val if pixel_val else 5)

    shore_point_layer.StartTransaction()
    for fid, risk_val_list in shore_point_feature_risk_map.items():
        min_rank = 5.0
        sum_sq_rank = 0.0
        for risk_val in risk_val_list:
            if risk_val < min_rank:
                min_rank = risk_val
            sum_sq_rank += (5 - risk_val)**2

        if sum_sq_rank > 0:
            r_hab_val = max(
                1.0, 4.8 - 0.5 * (
                    (1.5 * (5 - min_rank))**2 + sum_sq_rank -
                    (5 - min_rank)**2)**0.5)
        else:
            r_hab_val = 5.0

        shore_feature = shore_point_layer.GetFeature(fid)
        shore_feature.SetField(target_fieldname, float(r_hab_val))
        shore_point_layer.SetFeature(shore_feature)

    shore_point_layer.CommitTransaction()
    shore_point_layer = None
    shore_point_vector = None
    hab_effective_raster = None
    hab_effective_band = None
    retrying_rmtree(tmp_working_dir)


@retrying.retry(
    wait_exponential_multiplier=100, wait_exponential_max=2000,
    stop_max_attempt_number=5)
def retrying_rmtree(dir_path):
    """Remove `dir_path` but try a few times."""
    try:
        shutil.rmtree(dir_path)
    except Exception:
        LOGGER.exception('unable to remove %s' % dir_path)
        raise


def calculate_utm_srs(lng, lat):
    """Calculate UTM SRS from the lng/lat point given.

    Parameters:
        lng (float): longitude point.
        lat (float): latitude point.

    Returns:
        osr.SpatialReference in the UTM zone that contains the point (lng, lat)

    """
    utm_code = (math.floor((lng+180)/6) % 60) + 1
    lat_code = 6 if lat > 0 else 7
    epsg_code = int('32%d%02d' % (lat_code, utm_code))
    utm_srs = osr.SpatialReference()
    utm_srs.ImportFromEPSG(epsg_code)
    return utm_srs


def clip_geometry(
        bounding_box_coords, base_srs, target_srs, ogr_geometry_type,
        global_geom_strtree, target_vector_path):
    """Clip geometry in `global_geom_strtree` to bounding box.

    Parameters:
        bounding_box_coords (list): a list of bounding box coordinates in
            the same coordinate system as the geometry in
            `global_geom_strtree`.
        target_srs (osr.SpatialReference): target spatial reference for
            creating the target vector.
        ogr_geometry_type (ogr.wkb[TYPE]): geometry type to create for the
            target vector.
        global_geom_strtree (shapely.strtree.STRtree): an rtree loaded with
            geometry to query via bounding box. Each geometry will contain
            parameters `field_val_map` and `prep` that have values to copy to
            `target_fieldname` and used to quickly query geometry. Main object
            will have `field_name_type_list` field used to describe the
            original field name/types.
        target_vector_path (str): path to vector to create that will contain
            locally projected geometry clipped to the given bounding box.

    Returns:
        None.

    """
    gpkg_driver = ogr.GetDriverByName("GPKG")
    vector = gpkg_driver.CreateDataSource(
        target_vector_path)
    layer = vector.CreateLayer(
        os.path.splitext(os.path.basename(target_vector_path))[0],
        target_srs, ogr_geometry_type)
    for field_name, field_type in global_geom_strtree.field_name_type_list:
        layer.CreateField(ogr.FieldDefn(field_name, field_type))
    layer_defn = layer.GetLayerDefn()
    base_to_target_transform = osr.CoordinateTransformation(
        base_srs, target_srs)

    bounding_box = shapely.geometry.box(*bounding_box_coords)

    possible_geom_list = global_geom_strtree.query(bounding_box)
    LOGGER.debug('possible intersections %d', len(possible_geom_list))
    if not possible_geom_list:
        layer = None
        vector = None
        raise ValueError('no data intersects this box')
    for geom in possible_geom_list:
        clipped_line = bounding_box.intersection(geom)
        clipped_line_geom = ogr.CreateGeometryFromWkb(clipped_line.wkb)
        error_code = clipped_line_geom.Transform(base_to_target_transform)
        if error_code:
            raise RuntimeError(error_code)
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(clipped_line_geom.Clone())
        for field_name, _ in global_geom_strtree.field_name_type_list:
            feature.SetField(
                field_name, geom.field_val_map[field_name])
        layer.CreateFeature(feature)


def sample_line_to_points(
        line_vector_path, target_point_path, point_step_size):
    """Sample lines in line vector to points along the path.

    Parameters:
        line_vector_path (str): path to line based vector.
        target_point_path (str): created by this function. A GPKG that is in
            the same projection as `line_vector` where points lie on those line
            segments spaced no less than `point_step_size` apart.
        point_step_size (float): step size in projected units of `line_vector`
            for points to drop along the line segments.

    Returns:
        None.

    """
    line_vector = gdal.OpenEx(line_vector_path)
    line_layer = line_vector.GetLayer()
    layer_name = os.path.splitext(os.path.basename(line_vector_path))[0]

    gpkg_driver = ogr.GetDriverByName('GPKG')
    if os.path.exists(target_point_path):
        os.remove(target_point_path)
    point_vector = gpkg_driver.CreateDataSource(target_point_path)
    point_layer = point_vector.CreateLayer(
        layer_name, line_layer.GetSpatialRef(), ogr.wkbPoint,
        ['OVERWRITE=YES'])
    point_defn = point_layer.GetLayerDefn()
    for feature in line_layer:
        current_distance = 0.0
        line_geom = feature.GetGeometryRef()
        line = shapely.wkb.loads(line_geom.ExportToWkb())
        while current_distance < line.length:
            new_point = line.interpolate(current_distance)
            current_distance += point_step_size
            new_point_feature = ogr.Feature(point_defn)
            new_point_geom = ogr.CreateGeometryFromWkb(new_point.wkb)
            new_point_feature.SetGeometry(new_point_geom)
            point_layer.CreateFeature(new_point_feature)

    point_layer = None
    point_vector = None
    line_layer = None
    line_vector = None


def calculate_relief(
        shore_point_vector_path, dem_path, target_fieldname):
    """Calculate DEM relief as average coastal land area within 5km.

    Parameters:
        shore_point_vector_path (string):  path to a point shapefile to
            for relief point analysis.
        dem_path (string): path to a DEM raster projected in local coordinates.
        target_fieldname (string): this field name will be added to
            `shore_point_vector_path` and filled with Relief values.

    Returns:
        None.

    """
    try:
        shore_point_vector = gdal.OpenEx(
            shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
        shore_point_layer = shore_point_vector.GetLayer()
        relief_field = ogr.FieldDefn(target_fieldname, ogr.OFTReal)
        relief_field.SetPrecision(5)
        relief_field.SetWidth(24)
        shore_point_layer.CreateField(relief_field)
        dem_info = pygeoprocessing.get_raster_info(dem_path)

        tmp_working_dir = tempfile.mkdtemp(
            prefix='calculate_relief_',
            dir=os.path.dirname(shore_point_vector_path))

        dem_nodata = dem_info['nodata'][0]

        def zero_negative_values(depth_array):
            valid_mask = depth_array != dem_nodata
            result_array = numpy.empty_like(depth_array)
            result_array[:] = dem_nodata
            result_array[valid_mask] = 0
            result_array[depth_array > 0] = depth_array[depth_array > 0]
            return result_array

        positive_dem_path = os.path.join(
            tmp_working_dir, 'positive_dem.tif')

        pygeoprocessing.raster_calculator(
            [(dem_path, 1)], zero_negative_values,
            positive_dem_path, gdal.GDT_Int16, dem_nodata)

        # convolve over a 5km radius
        dem_pixel_size = dem_info['pixel_size']
        kernel_radius = (
            abs(RELIEF_SAMPLE_DISTANCE // dem_pixel_size[0]),
            abs(RELIEF_SAMPLE_DISTANCE // dem_pixel_size[1]))

        kernel_filepath = os.path.join(
            tmp_working_dir, 'averaging_kernel.tif')
        create_averaging_kernel_raster(kernel_radius, kernel_filepath)

        relief_path = os.path.join(tmp_working_dir, 'relief.tif')
        pygeoprocessing.convolve_2d(
            (positive_dem_path, 1), (kernel_filepath, 1), relief_path)
        relief_raster = gdal.Open(relief_path)
        relief_band = relief_raster.GetRasterBand(1)

        inv_gt = gdal.InvGeoTransform(dem_info['geotransform'])

        shore_point_layer.ResetReading()
        for point_feature in shore_point_layer:
            point_geometry = point_feature.GetGeometryRef()
            point_x, point_y = point_geometry.GetX(), point_geometry.GetY()
            point_geometry = None

            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(inv_gt, point_x, point_y)]

            try:
                pixel_value = relief_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception(
                    'relief_band size %d %d', relief_band.XSize,
                    relief_band.YSize)
                raise
            # Make relief "negative" so when we histogram it for risk a
            # "higher" value will show a lower risk.
            point_feature.SetField(target_fieldname, -float(pixel_value))
            shore_point_layer.SetFeature(point_feature)

        shore_point_layer.SyncToDisk()
        shore_point_layer = None
        shore_point_vector = None
        relief_raster = None
        relief_band = None

        try:
            retrying_rmtree(tmp_working_dir)
        except OSError:
            LOGGER.warning('unable to rm %s' % tmp_working_dir)

    except Exception:
        LOGGER.exception('error in relief calc')
        raise


def clip_and_reproject_raster(
        base_raster_path, target_raster_path, target_srs_wkt,
        target_bounding_box, edge_buffer, resample_method,
        reproject_bounding_box, target_pixel_size):
    """Clip and reproject base to target raster.

    Parameters:
        base_raster_path (str): path to the raster to clip from.
        target_raster_path (str): path to target raster that is a clip from
            base projected in `target_srs_wkt` coordinate system.
        target_srs_wkt (str): spatial reference of target coordinate system in
            wkt.
        target_bounding_box (list): List of float describing target bounding
            box in base coordinate system as [minx, miny, maxx, maxy].
        edge_buffer (float): amount to extend sides of bounding box in target
            coordinate system units.
        resample_method (str): one of
            "near|bilinear|cubic|cubicspline|lanczos|mode".
        reproject_bounding_box (bool): If true, project `target_bounding_box`
            from base coordinate system to `target_srs_wkt`.
        target_pixel_size (float): desired target pixel size in projected
            coordinates.

    Returns:
        None.

    """
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    bb_centroid = (
        (target_bounding_box[0]+target_bounding_box[2])/2,
        (target_bounding_box[1]+target_bounding_box[3])/2)

    if reproject_bounding_box:
        local_bounding_box = pygeoprocessing.transform_bounding_box(
            target_bounding_box, base_raster_info['projection'],
            target_srs_wkt, edge_samples=11)
    else:
        local_bounding_box = target_bounding_box
        base_srs = osr.SpatialReference()
        base_srs.ImportFromWkt(base_raster_info['projection'])
        target_srs = osr.SpatialReference()
        target_srs.ImportFromWkt(target_srs_wkt)
        target_to_base_transform = osr.CoordinateTransformation(
            target_srs, base_srs)
        point = ogr.CreateGeometryFromWkt("POINT (%f %f)" % bb_centroid)
        point.Transform(target_to_base_transform)
        bb_centroid = (point.GetX(), point.GetY())

    buffered_bounding_box = [
        local_bounding_box[0]-edge_buffer,
        local_bounding_box[1]-edge_buffer,
        local_bounding_box[2]+edge_buffer,
        local_bounding_box[3]+edge_buffer,
    ]

    # target_pixel_size = estimate_projected_pixel_size(
    #     base_raster_path, bb_centroid, target_srs_wkt)
    pygeoprocessing.warp_raster(
        base_raster_path, target_pixel_size, target_raster_path,
        resample_method, target_bb=buffered_bounding_box,
        target_sr_wkt=target_srs_wkt,
        working_dir=os.path.dirname(target_raster_path))


def clip_raster(
        base_raster_path, target_raster_path,
        target_bounding_box, edge_buffer):
    """Clip and reproject base to target raster.

    Parameters:
        base_raster_path (str): path to the raster to clip from.
        target_raster_path (str): path to target raster that is a clip from
            base projected in `target_srs_wkt` coordinate system.
        target_srs_wkt (str): spatial reference of target coordinate system in
            wkt.
        target_bounding_box (list): List of float describing target bounding
            box in base coordinate system as [minx, miny, maxx, maxy].
        edge_buffer (float): amount to extend sides of bounding box in target
            coordinate system units.
        resample_method (str): one of
            "near|bilinear|cubic|cubicspline|lanczos|mode".
        target_bounding_box (bool): If True, assumes bounding box is in
            base coordinate system and will transform it to target.

    Returns:
        None.

    """
    buffered_bounding_box = [
        target_bounding_box[0]-edge_buffer,
        target_bounding_box[1]-edge_buffer,
        target_bounding_box[2]+edge_buffer,
        target_bounding_box[3]+edge_buffer,
    ]

    base_raster = gdal.OpenEx(base_raster_path)
    gdal.Translate(
        target_raster_path, base_raster,
        format='GTiff',
        creationOptions=[
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
            'BLOCKXSIZE=256', 'BLOCKYSIZE=256'],
        outputBounds=buffered_bounding_box,
        callback=pygeoprocessing._make_logger_callback(
            "Translate %.1f%% complete"))
    base_raster = None


def estimate_projected_pixel_size(
        base_raster_path, sample_point, target_srs_wkt):
    """Estimate the pixel size of raster if projected in `target_srs_wkt`.

    Parameters:
        base_raster_path (str): path to a raster in some coordinate system.
        sample_point (list): [x, y] coordinate in base coordinate system of
            point to estimate projected pixel size around.
        target_srs_wkt (str): desired target coordinate system in wkt for
            estimate pixel size.

    Returns:
        None.

    """
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    base_pixel_size = base_raster_info['pixel_size']
    raster_center_pixel_bb = [
        sample_point[0] - abs(base_pixel_size[0]/2),
        sample_point[1] - abs(base_pixel_size[1]/2),
        sample_point[0] + abs(base_pixel_size[0]/2),
        sample_point[1] + abs(base_pixel_size[1]/2),
    ]
    pixel_bb = pygeoprocessing.transform_bounding_box(
        raster_center_pixel_bb, base_raster_info['projection'], target_srs_wkt)
    # x goes to the right, y goes down
    estimated_pixel_size = [
        pixel_bb[2]-pixel_bb[0],
        pixel_bb[1]-pixel_bb[3]]
    return estimated_pixel_size


def create_averaging_kernel_raster(radius_in_pixels, kernel_filepath):
    """Create a flat raster kernel with a 2d radius given.

    Parameters:
        radius_in_pixels (tuple): the (x/y) distance of the averaging kernel.
        kernel_filepath (string): The path to the file on disk where this
            kernel should be stored.  If this file exists, it will be
            overwritten.

    Returns:
        None

    """
    driver = gdal.GetDriverByName('GTiff')
    LOGGER.debug(radius_in_pixels)
    kernel_raster = driver.Create(
        kernel_filepath, int(2*radius_in_pixels[0]),
        int(2*radius_in_pixels[1]), 1, gdal.GDT_Float32)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([1, 0.1, 0, 1, 0, -0.1])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(-9999)

    n_cols = kernel_raster.RasterXSize
    n_rows = kernel_raster.RasterYSize
    iv, jv = numpy.meshgrid(range(n_rows), range(n_cols), indexing='ij')

    cx = n_cols / 2.0
    cy = n_rows / 2.0

    kernel_array = numpy.where(
        (cx-jv)**2 / radius_in_pixels[0] +
        (cy-iv)**2 / radius_in_pixels[1] <= 1.0, 1.0, 0.0)
    LOGGER.debug(kernel_array)

    # normalize
    kernel_array /= numpy.sum(kernel_array)
    kernel_band.WriteArray(kernel_array)
    kernel_band = None
    kernel_raster = None


def vector_to_lines(base_vector_path, target_line_vector_path):
    """Convert polygon vector to list of lines.

    Parameters:
        base_vector_path (str): path to polygon vector.
        target_line_vector_path (str): created by this file all polygons are
            converted to their line boundary equivalents.

    Returns:
        None.

    """
    # explode landmass into lines for easy intersection
    base_vector = gdal.OpenEx(base_vector_path, gdal.OF_VECTOR)
    base_layer = base_vector.GetLayer()

    gpkg_driver = ogr.GetDriverByName('GPKG')
    line_vector = gpkg_driver.CreateDataSource(
        target_line_vector_path)
    line_layer = line_vector.CreateLayer(
        target_line_vector_path, base_layer.GetSpatialRef(), ogr.wkbLineString)
    line_vector_defn = line_layer.GetLayerDefn()

    line_layer.StartTransaction()
    for base_feature in base_layer:
        base_shapely = shapely.wkb.loads(
            base_feature.GetGeometryRef().ExportToWkb())
        for line in geometry_to_lines(base_shapely):
            segment_feature = ogr.Feature(line_vector_defn)
            segement_geometry = ogr.Geometry(ogr.wkbLineString)
            segement_geometry.AddPoint(*line.coords[0])
            segement_geometry.AddPoint(*line.coords[1])
            segment_feature.SetGeometry(segement_geometry)
            line_layer.CreateFeature(segment_feature)
    line_layer.CommitTransaction()
    line_layer = None
    line_vector = None
    base_vector = None
    base_layer = None


def geometry_to_lines(geometry):
    """Convert a geometry object to a list of lines."""
    if geometry.type == 'Polygon':
        return polygon_to_lines(geometry)
    elif geometry.type == 'MultiPolygon':
        line_list = []
        for geom in geometry.geoms:
            line_list.extend(geometry_to_lines(geom))
        return line_list
    else:
        return []


def polygon_to_lines(geometry):
    """Return a list of shapely lines given higher order shapely geometry."""
    line_list = []
    last_point = geometry.exterior.coords[0]
    for point in geometry.exterior.coords[1::]:
        if point == last_point:
            continue
        line_list.append(shapely.geometry.LineString([last_point, point]))
        last_point = point
    line_list.append(shapely.geometry.LineString([
        last_point, geometry.exterior.coords[0]]))
    for interior in geometry.interiors:
        last_point = interior.coords[0]
        for point in interior.coords[1::]:
            if point == last_point:
                continue
            line_list.append(shapely.geometry.LineString([last_point, point]))
            last_point = point
        line_list.append(shapely.geometry.LineString([
            last_point, interior.coords[0]]))
    return line_list


def compute_wave_height(Un, Fn, dn):
    """Compute Wave Height by User Guide eq 10.

    This equation may not be suitable for wind speed values < 1 m/s
    The WWIII database tends to include some 0s, otherwise values > 2.

    Parameters:
        Un (float): wind velocity in meters per second.
        Fn (float): fetch ray length in meters.
        dn (float): water depth in negative meters.

    Returns:
        Float: Wave height in meters

    """
    if Un < 1.0:
        LOGGER.warning(
            'Found wind velocity of %.2f, '
            'using 1.0m/s in wave height calculation instead' % Un)
        Un = 1.0
    g = 9.81
    dn = -dn
    ds = g*dn/Un**2
    Fs = g*Fn/Un**2
    A = numpy.tanh(0.343*ds**1.14)
    B = numpy.tanh(4.41e-4*Fs**0.79/A)
    H_n = (0.24*Un**2/g)*(A*B)**0.572
    return H_n


def compute_wave_period(Un, Fn, dn):
    """Compute Wave Period by User Guide eq 10.

    This equation may not be suitable for wind speed values < 1 m/s
    The WWIII database tends to include some 0s, otherwise values > 2.

    Parameters:
        Un (float): wind velocity in meters per second.
        Fn (float): fetch ray length in meters.
        dn (float): water depth in negative meters.

    Returns:
        Float: Wave period in seconds

    """
    # This equation may not be suitable for wind speed values < 1 m/s
    # The WWIII database tends to include some 0s, otherwise values > 2
    if Un < 1.0:
        LOGGER.warning(
            'Found wind velocity of %.2f, '
            'using 1.0m/s in wave height calculation instead' % Un)
        Un = 1.0
    g = 9.81
    dn = -dn
    ds = g*dn/Un**2
    Fs = g*Fn/Un**2
    A = numpy.tanh(0.1*ds**2.01)
    B = numpy.tanh(2.77e-7*Fs**1.45/A)
    T_n = 7.69*Un/g*(A*B)**0.187
    return T_n


def merge_masks_op(mask_a, mask_b):
    return (mask_a == 1) & (mask_b == 1)


def merge_cv_points(cv_vector_queue, target_cv_vector_path):
    """Merge vectors in `cv_vector_queue` into single vector.

    Parameters:
        cv_vector_queue (multiprocessing.Processing): a queue containing
            paths to CV workspace point vectors. Terminated with
            `STOP_SENTINEL`.
        target_cv_vector_path (str): path to a point vector created by this
            function.

    Returns:
        None.

    """
    gpkg_driver = osr.GetDriverByName('GPKG')
    target_cv_vector = gpkg_driver.CreateDataSource(target_cv_vector_path)
    layer_name = os.path.basename(os.path.splitext(target_cv_vector_path)[0])
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    target_cv_layer = (
        target_cv_vector.CreateLayer(layer_name, wgs84_srs, ogr.wkbPoint))
    fields_to_copy = [
        'Rgeomorphology', 'surge', 'ew', 'rei', 'slr', 'Rhab', 'relief']
    for field_id in fields_to_copy:
        target_cv_layer.CreateField(ogr.FieldDefn(field_id, ogr.OFTReal))
    target_cv_layer_defn = target_cv_layer.GetLayerDefn()

    target_cv_layer.StartTransaction()
    while True:
        cv_vector_path = cv_vector_queue.get()
        if cv_vector_path == STOP_SENTINEL:
            break
        cv_vector = gdal.OpenEx(cv_vector_path, gdal.OF_VECTOR)
        cv_layer = cv_vector.GetLayer()
        cv_projection = cv_layer.GetSpatialRef()
        base_to_target_transform = osr.CoordinateTransformation(
            cv_projection, wgs84_srs)

        for cv_feature in cv_layer:
            cv_geom = cv_feature.GetGeometryRef().Clone()
            _ = cv_geom.Transform(base_to_target_transform)
            target_feature = ogr.Feature(target_cv_layer_defn)
            target_feature.SetGeometry(cv_geom)
            for field_id in fields_to_copy:
                target_feature.SetField(cv_feature.GetField(field_id))
            target_cv_layer.CreateFeature(target_feature)
        cv_feature = None
        cv_geom = None
        cv_layer = None
        cv_vector = None
    target_cv_layer.CommitTransaction()


if __name__ == '__main__':
    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR, GRID_WORKSPACE_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(CHURN_DIR, -1, 5.0)

    for zip_url in [LS_POPULATION_RASTER_URL]:
        target_token_path = os.path.join(
            CHURN_DIR, os.path.basename(os.path.splitext(zip_url)[0]))
        download_and_unzip_task = task_graph.add_task(
            func=download_and_unzip,
            args=(zip_url, ECOSHARD_DIR, target_token_path),
            target_path_list=[target_token_path],
            task_name='download and unzip %s' % zip_url)

    local_data_path_map = {}

    for data_id, ecoshard_url in GLOBAL_DATA_URL_MAP.items():
        local_ecoshard_path = os.path.join(
            ECOSHARD_DIR, os.path.basename(ecoshard_url))
        download_task = task_graph.add_task(
            func=ecoshard.download_url,
            args=(ecoshard_url, local_ecoshard_path),
            target_path_list=[local_ecoshard_path],
            task_name='download %s' % local_ecoshard_path)
        local_data_path_map[data_id] = local_ecoshard_path

    ls_population_raster_path = os.path.join(ECOSHARD_DIR, 'lspop2017')

    global_wwiii_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(
            os.path.splitext(GLOBAL_WWIII_GZ_URL)[0]))
    download_and_ungzip_global_wwiii_task = task_graph.add_task(
        func=download_and_ungzip,
        args=(GLOBAL_WWIII_GZ_URL, global_wwiii_vector_path),
        target_path_list=[global_wwiii_vector_path],
        task_name='download %s' % global_wwiii_vector_path)

    shore_buffer_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(BUFFER_VECTOR_URL))
    download_buffer_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(BUFFER_VECTOR_URL, shore_buffer_vector_path),
        target_path_list=[shore_buffer_vector_path],
        task_name='download global_vector')

    task_graph.join()  # wait for everything to download

    lulc_shore_mask_raster_path = os.path.join(
        CHURN_DIR, 'lulc_masked_by_shore.tif')
    mask_lulc_by_shore_task = task_graph.add_task(
        func=pygeoprocessing.mask_raster,
        args=((local_data_path_map['lulc'], 1), shore_buffer_vector_path,
              lulc_shore_mask_raster_path),
        target_path_list=[lulc_shore_mask_raster_path],
        ignore_path_list=[shore_buffer_vector_path],  # ignore mod by opening
        dependent_task_list=[download_buffer_task],
        task_name='mask shore')

    # each value in `risk_distance_to_lulc_code` can be lumped into one
    # type.
    risk_distance_to_lulc_code = collections.defaultdict(list)
    for lulc_code, risk_distance in LULC_CODE_TO_HAB_MAP.items():
        risk_distance_to_lulc_code[risk_distance].append(lulc_code)

    # this maps all the same type of codes together
    lulc_code_to_reclass_value = {}
    habitat_raster_risk_map = dict(HABITAT_VECTOR_PATH_MAP)
    for risk_distance_tuple, lulc_code_list in sorted(
            risk_distance_to_lulc_code.items()):
        if risk_distance_tuple[0] == 0:
            LOGGER.info('skipping hab tuple %s', str(risk_distance_tuple))
            continue
        # reclassify landcover map to be ones everywhere for `lulc_code_list`
        reclass_map = {}
        for lulc_code in LULC_CODE_TO_HAB_MAP:
            if lulc_code in lulc_code_list:
                reclass_map[lulc_code] = 1
            else:
                reclass_map[lulc_code] = 0

        risk_distance_mask_path = os.path.join(
            CHURN_DIR, '%s_%s_mask.tif' % risk_distance_tuple)

        risk_distance_mask_task = task_graph.add_task(
            func=pygeoprocessing.reclassify_raster,
            args=(
                (lulc_shore_mask_raster_path, 1), reclass_map,
                risk_distance_mask_path, gdal.GDT_Byte, 0),
            target_path_list=[risk_distance_mask_path],
            dependent_task_list=[mask_lulc_by_shore_task],
            task_name='map distance types %s' % str(risk_distance_tuple))
        habitat_raster_risk_map[risk_distance_tuple] = (
            risk_distance_mask_path, risk_distance_tuple[0],
            risk_distance_tuple[1])

    # TODO: merge mangroves and onshore forest
    # 'mangroves'
    # (1, 2000) -- forest risk/dist
    merged_forest_mangrove_raster_path = os.path.join(
        CHURN_DIR, 'merged_forest_mangrove.tif')
    aligned_raster_path_list = [
        os.path.join(CHURN_DIR, x) for x in ['a.tif', 'b.tif']]
    habitat_raster_info = pygeoprocessing.get_raster_info(
        habitat_raster_risk_map[(1, 2000)][0])
    align_task = task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            [habitat_raster_risk_map[(1, 2000)][0],
             local_data_path_map['mangroves']], aligned_raster_path_list,
            ['near', 'near'],
            habitat_raster_info['pixel_size'], 'union'),
        target_path_list=aligned_raster_path_list,
        task_name='align mangrove and forest')

    mask_task = task_graph.add_task(
        func=pygeoprocessing.raster_calculator,
        args=(
            ((aligned_raster_path_list[0], 1),
             (aligned_raster_path_list[1], 1)), merge_masks_op,
            merged_forest_mangrove_raster_path, 0, gdal.GDT_Byte),
        target_path_list=[merged_forest_mangrove_raster_path],
        dependent_task_list=[align_task],
        task_name='merge mangrove/forest masks')

    task_graph.join()
    task_graph.close()

    del habitat_raster_risk_map[(1, 2000)]
    habitat_raster_risk_map['mangroves'] = (
        merged_forest_mangrove_raster_path, 1, 1000)

    shore_grid_vector = gdal.OpenEx(
        local_data_path_map['shore_grid'], gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    bb_work_queue = multiprocessing.Queue()
    cv_point_complete_queue = multiprocessing.Queue()

    cv_grid_worker_thread = threading.Thread(
        target=cv_grid_worker,
        args=(
            bb_work_queue,
            cv_point_complete_queue,
            local_data_path_map['landmass'],
            local_data_path_map['geomorphology'],
            local_data_path_map['slr'],
            local_data_path_map['dem'],
            global_wwiii_vector_path,
            habitat_raster_risk_map,
            ))

    for path in [
            ls_population_raster_path,
            local_data_path_map['lulc'],
            global_wwiii_vector_path,
            local_data_path_map['landmass'],
            local_data_path_map['shore_grid']]:
        LOGGER.info('%s: %s' % (os.path.exists(path), path))

    shore_grid_vector = gdal.OpenEx(
        local_data_path_map['shore_grid'], gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    for index, shore_grid_feature in enumerate(shore_grid_layer):
        shore_grid_geom = shore_grid_feature.GetGeometryRef()
        boundary_box = shapely.wkb.loads(shore_grid_geom.ExportToWkb())
        LOGGER.debug(boundary_box.bounds)
        bb_work_queue.put((index, boundary_box.bounds))
        if index > 10:
            break

    bb_work_queue.put(STOP_SENTINEL)
    cv_grid_worker_thread.start()

    merge_cv_points_thread = threading.Thread(
        target=merge_cv_points,
        args=(cv_point_complete_queue, TARGET_CV_VECTOR_PATH))
    merge_cv_points_thread.start()

    cv_grid_worker_thread.join()
    # when workers are complete signal merger complete
    cv_point_complete_queue.put(STOP_SENTINEL)

    merge_cv_points_thread.join()

    LOGGER.debug('cv grid joined')
