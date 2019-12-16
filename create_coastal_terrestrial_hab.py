"""Buffer CV habitat."""
import sys
import logging
import os
import urllib.request

from osgeo import gdal
import pygeoprocessing
import taskgraph

LULC_RASTER_URL = 'https://storage.googleapis.com/ipbes-ndr-ecoshard-data/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'
BUFFER_VECTOR_URL = 'https://storage.googleapis.com/ecoshard-root/working_shards/buffered_global_shore_5km_md5_a68e1049c1c03673add014cd29b7b368.gpkg'


# I got this from https://docs.google.com/spreadsheets/d/1pYNWwPBqYYZ4tdJC3zaAZ8Z-CMOU2bZnB-ZvhKzDQlU/edit#gid=0

LULC_CODE_TO_HAB_MAP = {
    50: 1, 60: 1, 61: 1, 62: 1, 70: 1, 71: 1, 72: 1, 80: 1, 81: 1, 82: 1,
    90: 1, 100: 1, 110: 2, 120: 2, 121: 2, 122: 2, 130: 2, 140: 2, 150: 3,
    152: 4, 153: 4, 160: 5, 170: 5, 180: 5, 190: 5, 0: 0, 10: 0, 11: 0, 12: 0,
    20: 0, 30: 0, 40: 0, 200: 0, 201: 0, 202: 0, 210: 0, 220: 0,
    151: 3}

WORKSPACE_DIR = 'buffer_habitat_workspace'
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
LOGGING_LEVEL = logging.DEBUG

logging.basicConfig(
    level=LOGGING_LEVEL,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def main():
    """Entry point."""
    for dir_path in [ECOSHARD_DIR, CHURN_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, -1, 5)

    buffer_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(BUFFER_VECTOR_URL))
    download_buffer_task = task_graph.add_task(
        func=urllib.request.urlretrieve,
        args=(BUFFER_VECTOR_URL, buffer_vector_path),
        target_path_list=[buffer_vector_path],
        task_name='download global_vector')

    lulc_raster_path = os.path.join(
        WORKSPACE_DIR, os.path.basename(LULC_RASTER_URL))
    download_lulc_task = task_graph.add_task(
        func=urllib.request.urlretrieve,
        args=(LULC_RASTER_URL, lulc_raster_path),
        target_path_list=[lulc_raster_path],
        task_name='download lulc raster')

    shore_mask_path = os.path.join(CHURN_DIR, 'shore_mask.tif')
    mask_shore_task = task_graph.add_task(
        func=make_mask,
        args=(lulc_raster_path, buffer_vector_path, shore_mask_path),
        target_path_list=[shore_mask_path],
        dependent_task_list=[download_buffer_task, download_lulc_task],
        task_name='mask shore')

    cv_habitat_map_raster_path = os.path.join(CHURN_DIR, 'cv_habitat_map.tif')
    lulc_nodata = 0
    task_graph.add_task(
        func=pygeoprocessing.reclassify_raster,
        args=(
            (lulc_raster_path, 1), LULC_CODE_TO_HAB_MAP,
            cv_habitat_map_raster_path, gdal.GDT_Byte, lulc_nodata),
        target_path_list=[cv_habitat_map_raster_path],
        dependent_task_list=[download_lulc_task],
        task_name='map lulc to cv types')

    task_graph.join()
    task_graph.close()


def make_mask(base_raster_path, mask_vector_path, target_raster_path):
    """Mask vector onto target using base as the reference."""
    pygeoprocessing.new_raster_from_base(
        base_raster_path, target_raster_path, gdal.GDT_Byte, [0],
        fill_value_list=[0])
    pygeoprocessing.rasterize(
        mask_vector_path, target_raster_path, burn_values=[1])


if __name__ == '__main__':
    main()
