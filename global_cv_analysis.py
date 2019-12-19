"""Global CV Analysis for CNC.

Design doc is here:

https://docs.google.com/document/d/18AcJM-rXeIYgEsmqlaUwdtm7gdWLiaD6kkRkpARILlw/edit#heading=h.bbujb61ete53

"""
import collections
import gzip
import logging
import os
import sys
import zipfile

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pygeoprocessing
import shapely.geometry
import shapely.wkb
import shapely.strtree
import taskgraph

WORKING_DIR = "cnc_cv_workspace"
ECOSHARD_DIR = os.path.join(WORKING_DIR, 'ecoshard')
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
GLOBAL_WWIII_GZ_URL = (
    ECOSHARD_BUCKET_URL +
    'wave_watch_iii_md5_c8bb1ce4739e0a27ee608303c217ab5b.gpkg.gz')
GLOBAL_DEM_URL = (
    ECOSHARD_BUCKET_URL +
    'global_dem_md5_22c5c09ac4c4c722c844ab331b34996c.tif')
LS_POPULATION_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')

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

WORKSPACE_DIR = 'global_cv_workspace'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')

SHORE_GRID_VECTOR_PATH = os.path.join(CHURN_DIR, 'shore_grid.gpkg')


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


if __name__ == '__main__':
    for dir_path in [WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR]:
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

    ls_population_raster_path = os.path.join(ECOSHARD_DIR, 'lspop2017')

    global_polygon_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(GLOBAL_POLYGON_URL))
    download_global_polygon_path = task_graph.add_task(
        func=ecoshard.download_url,
        args=(GLOBAL_POLYGON_URL, global_polygon_path),
        target_path_list=[global_polygon_path],
        task_name='download %s' % global_polygon_path)

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
        dependent_task_list=[download_buffer_task, download_lulc_task],
        task_name='mask shore')

    # each value in `risk_distance_to_lulc_code` can be lumped into one
    # type.
    risk_distance_to_lulc_code = collections.defaultdict(list)
    for lulc_code, risk_distance in LULC_CODE_TO_HAB_MAP.items():
        risk_distance_to_lulc_code[risk_distance].append(lulc_code)

    # this maps all the same type of codes together
    lulc_code_to_reclass_value = {}
    risk_distance_path_map = {}
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

    task_graph.join()
    task_graph.close()

    GLOBAL_AOI_WGS84_BB = [-179, -65, 180, 77]
    land_geometry_list = []

    global_polygon_vector = gdal.OpenEx(global_polygon_path, gdal.OF_VECTOR)
    global_polygon_layer = global_polygon_vector.GetLayer()
    feature_count = global_polygon_layer.GetFeatureCount()
    for index, global_polygon_feature in enumerate(global_polygon_layer):
        if index % 1000 == 0:
            LOGGER.debug(
                'adding %d of %d %.2f', index, feature_count,
                100.0 * index / feature_count)
        global_polygon_geom = global_polygon_feature.GetGeometryRef()
        global_polygon_shapely = (
            shapely.wkb.loads(global_polygon_geom.ExportToWkb()))
        global_polygon_prep = shapely.prepared.prep(global_polygon_shapely)
        global_polygon_shapely.prep = global_polygon_prep
        land_geometry_list.append(global_polygon_shapely)

    geometry_r_tree = shapely.strtree.STRtree(land_geometry_list)

    gpkg_driver = ogr.GetDriverByName('GPKG')
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    shore_grid_vector = gpkg_driver.CreateDataSource(
        SHORE_GRID_VECTOR_PATH)
    shore_grid_layer = shore_grid_vector.CreateLayer(
        os.path.splitext(os.path.basename(
            SHORE_GRID_VECTOR_PATH))[0], wgs84_srs, ogr.wkbPolygon)
    shore_grid_layer_defn = shore_grid_layer.GetLayerDefn()

    for lng in range(GLOBAL_AOI_WGS84_BB[0], GLOBAL_AOI_WGS84_BB[2]):
        for lat in range(GLOBAL_AOI_WGS84_BB[1], GLOBAL_AOI_WGS84_BB[3]):
            sample_box = shapely.geometry.box(lng, lat, lng+1, lat+1)
            intersection_list = geometry_r_tree.query(sample_box)
            for intersection_geom in intersection_list:
                if (intersection_geom.prep.intersects(sample_box) and
                        not intersection_geom.prep.contains(sample_box)):
                    LOGGER.debug('edge box: %s', str(sample_box))
                    box_feature = ogr.Feature(shore_grid_layer_defn)
                    sample_box_geom = ogr.CreateGeometryFromWkb(sample_box.wkb)
                    box_feature.SetGeometry(sample_box_geom)
                    shore_grid_layer.CreateFeature(box_feature)

    for path in [
            ls_population_raster_path,
            geomorphology_vector_path,
            reefs_vector_path,
            mangroves_vector_path,
            seagrass_vector_path,
            global_saltmarsh_vector_path,
            lulc_raster_path,
            global_wwiii_vector_path,
            global_polygon_path]:
        LOGGER.info('%s: %s' % (os.path.exists(path), path))
