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
import sys
import threading
import zipfile

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pygeoprocessing
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

#[minx, miny, maxx, maxy].
GLOBAL_AOI_WGS84_BB = [-179, -65, 180, 77]

ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')
GLOBAL_POLYGON_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_global_polygon_simplified_geometries_'
    'md5_653118dde775057e24de52542b01eaee.gpkg')
GLOBAL_GEOMORPHOLOGY_ZIP_URL = (
    ECOSHARD_BUCKET_URL + 'SedType_md5_d670148183cc7f817e26570e77ca3449.zip')
GLOBAL_REEFS_ZIP_URL = (
    ECOSHARD_BUCKET_URL + 'Reefs_md5_50b159fff92762ffa5efa0c611b436ce.zip')
GLOBAL_MANGROVES_ZIP_URL = (
    ECOSHARD_BUCKET_URL +
    'GMW_001_GlobalMangroveWatch_2016-20191216T181028Z-001_'
    'md5_dd4b04fb9cac2314af99beed2b2c856a.zip')
GLOBAL_SEAGRASS_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_seagrass_valid_md5_e206dde7cc9b95ba9846efa12b63d333.gpkg')
GLOBAL_SALTMARSH_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_saltmarsh_valid_md5_56364edc15ab96d79b9fa08b12ec56ab.gpkg')

LULC_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_'
    'md5_1254d25f937e6d9bdee5779d377c5aa4.tif')
BUFFER_VECTOR_URL = (
    ECOSHARD_BUCKET_URL +
    'buffered_global_shore_5km_md5_a68e1049c1c03673add014cd29b7b368.gpkg')
SHORE_GRIDS_URL = (
    ECOSHARD_BUCKET_URL +
    'shore_grid_md5_1a09b69c0c548d3b25d7e14c3ddb60c9.gpkg')

GLOBAL_WWIII_GZ_URL = (
    ECOSHARD_BUCKET_URL +
    'wave_watch_iii_md5_c8bb1ce4739e0a27ee608303c217ab5b.gpkg.gz')
GLOBAL_DEM_URL = (
    ECOSHARD_BUCKET_URL +
    'global_dem_md5_22c5c09ac4c4c722c844ab331b34996c.tif')
LS_POPULATION_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')
SLR_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'MSL_Map_MERGED_Global_AVISO_NoGIA_Adjust_'
    'md5_3072845759841d0b2523d00fe9518fee.tif')

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
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ',
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

STOP_SENTINEL = 'STOP'

HABITAT_VECTOR_PATH_MAP = {
    'reefs': (
        os.path.join(ECOSHARD_DIR, 'Reefs', 'reef_500_poly.shp'), 1, 2000.0),
    'mangroves': (
        os.path.join(
            ECOSHARD_DIR, 'GMW_001_GlobalMangroveWatch_2016',
            '01_Data', 'GMW_2016_v2.shp'), 1, 1000.0),
    'saltmarsh': (
        os.path.join(
            ECOSHARD_DIR,
            'ipbes-cv_saltmarsh_valid_'
            'md5_56364edc15ab96d79b9fa08b12ec56ab.gpkg'), 2, 1000.0),
    'seagrass': (
        os.path.join(
            ECOSHARD_DIR,
            'ipbes-cv_seagrass_valid_'
            'md5_e206dde7cc9b95ba9846efa12b63d333.gpkg'), 4, 500.0)}


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


def build_rtree(vector_path, field_to_copy=None):
    """Build an rtree that generates geom and preped geometry.

    Parameters:
        vector_path (str): path to vector of geometry to build into
            r tree.
        field_to_copy (str): if not None, the value of this field per feature
            will be copied into the .field_val parameter of each object in the
            r-tree.

    Returns:
        strtree.STRtree object that will return shapely geometry objects
            with a .prep field that is prepared geomtry for fast testing and
            a .geom field that is the base gdal geometry.

    """
    geometry_prep_list = []
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    for feature in layer:
        feature_geom = feature.GetGeometryRef().Clone()
        feature_geom_shapely = shapely.wkb.loads(feature_geom.ExportToWkb())
        feature_geom_shapely.prep = shapely.prepared.prep(feature_geom_shapely)
        feature_geom_shapely.geom = feature_geom
        if field_to_copy:
            feature_geom_shapely.field_val = feature.GetField(field_to_copy)
        geometry_prep_list.append(feature_geom_shapely)
    r_tree = shapely.strtree.STRtree(geometry_prep_list)
    return r_tree


# Rgeomorphology
# Rrelief
# Rhab
# Rslr
# Rwind
# Rwave
# Rsurge

def cv_grid_worker(
        bb_work_queue,
        geomorphology_vector_path,
        global_dem_raster_path,
        habitat_vector_path_map,
        habitat_raster_path_map,
        slr_vector_path,
        ww_iii_vector_path,
        ):
    """Worker process to calculate CV for a grid.

    Parameters:
        bb_work_queue (multiprocessing.Queue): contains
            [minx, miny, maxx, maxy] bounding box values to be processed or
            `STOP_SENTINEL` values to indicate the worker should be terminated.
        habitat_vector_path_map (dict): maps a habitat id to a
            (path, risk, dist(m)) tuple. These habitats should be used in the
            calculation of Rhab.
        habitat_raster_path_map (dict): mapt a habitat id to a
            (raster path, risk, dist(m)) tuple. These are the raster versions
            of habitats to use in Rhab.

    Returns:
        None.

    """
    geomorphology_rtree = build_rtree(geomorphology_vector_path, 'SEDTYPE')
    geomorphology_proj_wkt = pygeoprocessing.get_vector_info(
        geomorphology_vector_path)['projection']
    gegeomorphology_proj = osr.SpatialReference()
    gegeomorphology_proj.ImportFromWkt(geomorphology_proj_wkt)
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
        lng_min, lat_min, lng_max, lat_max = payload
        bounding_box = shapely.geometry.box(*payload)
        # create workspace
        workspace_dir = os.path.join(
            GRID_WORKSPACE_DIR, '%s_%s_%s_%s' % payload)

        try:
            os.makedirs(workspace_dir)
        except OSError:
            pass

        # task_graph = taskgraph.TaskGraph(workspace_dir, -1)

        utm_srs = calculate_utm_srs((lng_min+lng_max)/2, (lat_min+lat_max)/2)

        ### geomorphology
        local_geomorphology_vector_path = os.path.join(
            workspace_dir, 'geomorphology.gpkg')
        clip_geometry(
            bounding_box, utm_srs,
            ogr.wkbLineString, geomorphology_rtree, 'risk',
            local_geomorphology_vector_path)

        # Rrelief
        # Rhab
        # Rslr
        # Rwind
        # Rwave
        # Rsurge


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
        bounding_box, target_srs, ogr_geometry_type,
        global_geom_rtree, target_field_name, target_vector_path):
    """Clip geometry in `global_geom_rtree` to bounding box.

    Parameters:
        boudnging_box (shapely.Box): a shapely box in the same coordinate
            system as the geometry in `global_geom_rtree`.
        target_srs (osr.SpatialReference): target spatial reference for
            creating the target vector.
        ogr_geometry_type (ogr.wkb[TYPE]): geometry type to create for the
            target vector.
        global_geom_rtree (shapely.str_rtree): an rtree loaded with geometry
            to query via bounding box. Each geometry will contain parameters
            `field_val` and `prep` that have values to copy to
            `target_field_name` and used to quickly query geometry.
        target_field_name (str): field name to create in the target feature
            that will contain the value of .field_val
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
    layer.CreateField(
        ogr.FieldDefn(target_field_name, ogr.OFTReal))
    layer_defn = layer.GetLayerDefn()

    possible_geom_list = global_geom_rtree.query(bounding_box)
    if not possible_geom_list:
        raise ValueError('no data intersects this box')
    for geom in possible_geom_list:
        clipped_line = bounding_box.intersection(geom)
        clipped_line_geom = ogr.CreateGeometryFromWkb(clipped_line.wkb)
        error_code = clipped_line_geom.TransformTo(target_srs)
        if error_code:
            raise RuntimeError(error_code)
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(clipped_line_geom.Clone())
        feature.SetField('risk', geom.field_val)
        layer.CreateFeature(feature)


if __name__ == '__main__':
    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR, GRID_WORKSPACE_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(CHURN_DIR, -1, 5.0)

    for zip_url in [
            GLOBAL_GEOMORPHOLOGY_ZIP_URL, GLOBAL_REEFS_ZIP_URL,
            GLOBAL_MANGROVES_ZIP_URL, LS_POPULATION_URL]:
        target_token_path = os.path.join(
            CHURN_DIR, os.path.basename(os.path.splitext(zip_url)[0]))
        download_and_unzip_task = task_graph.add_task(
            func=download_and_unzip,
            args=(zip_url, ECOSHARD_DIR, target_token_path),
            target_path_list=[target_token_path],
            task_name='download and unzip %s' % zip_url)

    slr_raster_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(SLR_RASTER_URL))
    download_global_polygon_path = task_graph.add_task(
        func=ecoshard.download_url,
        args=(SLR_RASTER_URL, slr_raster_path),
        target_path_list=[slr_raster_path],
        task_name='download %s' % slr_raster_path)

    ls_population_raster_path = os.path.join(ECOSHARD_DIR, 'lspop2017')

    global_polygon_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(GLOBAL_POLYGON_URL))
    download_global_polygon_path = task_graph.add_task(
        func=ecoshard.download_url,
        args=(GLOBAL_POLYGON_URL, global_polygon_path),
        target_path_list=[global_polygon_path],
        task_name='download %s' % global_polygon_path)

    shore_grid_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(SHORE_GRIDS_URL))

    download_global_polygon_path = task_graph.add_task(
        func=ecoshard.download_url,
        args=(SHORE_GRIDS_URL, shore_grid_vector_path),
        target_path_list=[shore_grid_vector_path],
        task_name='download %s' % shore_grid_vector_path)

    global_dem_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(GLOBAL_DEM_URL))
    download_global_dem_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(GLOBAL_DEM_URL, global_dem_path),
        target_path_list=[global_dem_path],
        task_name='download %s' % global_dem_path)

    global_wwiii_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(
            os.path.splitext(GLOBAL_WWIII_GZ_URL)[0]))
    download_and_ungzip_global_wwiii_task = task_graph.add_task(
        func=download_and_ungzip,
        args=(GLOBAL_WWIII_GZ_URL, global_wwiii_vector_path),
        target_path_list=[global_wwiii_vector_path],
        task_name='download %s' % global_wwiii_vector_path)

    seagrass_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(GLOBAL_SEAGRASS_URL))
    download_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(GLOBAL_SEAGRASS_URL, seagrass_vector_path),
        target_path_list=[seagrass_vector_path],
        task_name='download seagrass task')

    global_saltmarsh_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(GLOBAL_SALTMARSH_URL))
    download_global_saltmarsh_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(GLOBAL_SALTMARSH_URL, global_saltmarsh_vector_path),
        target_path_list=[global_saltmarsh_vector_path],
        task_name='download global saltmarsh')

    mangroves_vector_path = os.path.join(
        ECOSHARD_DIR, 'GMW_001_GlobalMangroveWatch_2016', '01_Data',
        'GMW_2016_v2.shp')
    reefs_vector_path = os.path.join(
        ECOSHARD_DIR, 'Reefs', 'reef_500_poly.shp')
    geomorphology_vector_path = os.path.join(ECOSHARD_DIR, 'Sediments.shp')
    lulc_raster_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(LULC_RASTER_URL))

    download_lulc_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(LULC_RASTER_URL, lulc_raster_path),
        target_path_list=[lulc_raster_path],
        task_name='download lulc raster')

    shore_buffer_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(BUFFER_VECTOR_URL))
    download_buffer_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(BUFFER_VECTOR_URL, shore_buffer_vector_path),
        target_path_list=[shore_buffer_vector_path],
        task_name='download global_vector')

    lulc_shore_mask_raster_path = os.path.join(
        CHURN_DIR, 'lulc_masked_by_shore.tif')
    mask_lulc_by_shore_task = task_graph.add_task(
        func=pygeoprocessing.mask_raster,
        args=((lulc_raster_path, 1), shore_buffer_vector_path,
              lulc_shore_mask_raster_path),
        target_path_list=[lulc_shore_mask_raster_path],
        ignore_path_list=[shore_buffer_vector_path],  # ignore mod by opening
        dependent_task_list=[download_buffer_task, download_lulc_task],
        task_name='mask shore')

    # each value in `risk_distance_to_lulc_code` can be lumped into one
    # type.
    risk_distance_to_lulc_code = collections.defaultdict(list)
    for lulc_code, risk_distance in LULC_CODE_TO_HAB_MAP.items():
        risk_distance_to_lulc_code[risk_distance].append(lulc_code)

    # this maps all the same type of codes together
    lulc_code_to_reclass_value = {}
    habitat_raster_risk_map = {}
    for risk_distance_tuple, lulc_code_list in sorted(
            risk_distance_to_lulc_code.items()):
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
        habitat_raster_risk_map[risk_distance_tuple] = risk_distance_mask_path

    task_graph.join()
    task_graph.close()

    shore_grid_vector = gdal.OpenEx(shore_grid_vector_path, gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    bb_work_queue = multiprocessing.Queue()

    cv_grid_worker_thread = threading.Thread(
        target=cv_grid_worker,
        args=(
            bb_work_queue,
            geomorphology_vector_path,
            global_dem_path,
            HABITAT_VECTOR_PATH_MAP,
            habitat_raster_risk_map,
            slr_raster_path,
            global_wwiii_vector_path,
            ))

    for path in [
            ls_population_raster_path,
            geomorphology_vector_path,
            reefs_vector_path,
            mangroves_vector_path,
            seagrass_vector_path,
            global_saltmarsh_vector_path,
            lulc_raster_path,
            global_wwiii_vector_path,
            global_polygon_path,
            shore_grid_vector_path]:
        LOGGER.info('%s: %s' % (os.path.exists(path), path))

    shore_grid_vector = gdal.OpenEx(shore_grid_vector_path, gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    for shore_grid_feature in shore_grid_layer:
        shore_grid_geom = shore_grid_feature.GetGeometryRef()
        boundary_box = shapely.wkb.loads(shore_grid_geom.ExportToWkb())
        LOGGER.debug(boundary_box.bounds)
        bb_work_queue.put(boundary_box.bounds)

    bb_work_queue.put(STOP_SENTINEL)
    cv_grid_worker_thread.start()
    cv_grid_worker_thread.join()
    LOGGER.debug('cv grid joined')
