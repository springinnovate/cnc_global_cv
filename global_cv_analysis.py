"""Global CV Analysis for CNC.

Design doc is here:

https://docs.google.com/document/d/18AcJM-rXeIYgEsmqlaUwdtm7gdWLiaD6kkRkpARILlw/edit#heading=h.bbujb61ete53

"""
import logging
import os
import sys

from osgeo import gdal

WORKING_DIR = "cnc_cv_workspace"
ECOSHARD_DIR = os.path.join(WORKING_DIR, 'ecoshard')
_TARGET_NODATA = -1

_ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')

_GLOBAL_GEOMORPHOLOGY_URL = (
    _ECOSHARD_BUCKET_URL + 'SedType_md5_d670148183cc7f817e26570e77ca3449.zip')
_GLOBAL_REEFS_URL = (
    _ECOSHARD_BUCKET_URL + 'Reefs_md5_50b159fff92762ffa5efa0c611b436ce.zip')
_GLOBAL_MANGROVES_URL = (
    _ECOSHARD_BUCKET_URL +
    'GMW_001_GlobalMangroveWatch_2016-20191216T181028Z-001_'
    'md5_dd4b04fb9cac2314af99beed2b2c856a.zip')
_GLOBAL_SEAGRASS_URL = (
    _ECOSHARD_BUCKET_URL +
    'ipbes-cv_seagrass_valid_md5_e206dde7cc9b95ba9846efa12b63d333.gpkg')


# tuple form is (path, divide by area?, area to search, extra pixels to add)
_AGGREGATION_LAYER_MAP = {

    'pdn_gpw': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/gpw-v4-population-count-2015/gpw-v4-population-count_2015.tif"), True, None, 1e3, 0),
    'pdn_ssp1': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/Spatial_population_scenarios_GeoTIFF/SSP1_GeoTIFF/total/GeoTIFF/ssp1_2050.tif"), True, None, 1e3, 0),
    'pdn_ssp3': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/Spatial_population_scenarios_GeoTIFF/SSP3_GeoTIFF/total/GeoTIFF/ssp3_2050.tif"), True, None, 1e3, 0),
    'pdn_ssp5': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/Spatial_population_scenarios_GeoTIFF/SSP5_GeoTIFF/total/GeoTIFF/ssp5_2050.tif"), True, None, 1e3, 0),
    'pdn_2010': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/Spatial_population_scenarios_GeoTIFF/SSP1_GeoTIFF/total/GeoTIFF/ssp1_2010.tif"), True, None, 1e3, 0),
    '14bt_pop': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/gpw_v4_e_a000_014bt_2010_cntm_30_sec.tif"), True, None, 1e3, 0),
    '65plus_pop': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/gpw_v4_e_a065plusbt_2010_cntm_30_sec.tif"), True, None, 1e3, 0),
    'urbp_2015': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/GLOBIO4_landuse_10sec_tifs_20171207_Idiv/Current2015/Globio4_landuse_10sec_2015_cropint.tif"), False, [1, 190], 5e3, 0),
    'urbp_ssp1': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/GLOBIO4_landuse_10sec_tifs_20171207_Idiv/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif"), False, [1, 190], 5e3, 0),
    'urbp_ssp3': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/GLOBIO4_landuse_10sec_tifs_20171207_Idiv/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif"), False, [1, 190], 5e3, 0),
    'urbp_ssp5': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/GLOBIO4_landuse_10sec_tifs_20171207_Idiv/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif"), False, [1, 190], 5e3, 0),
    'SLRrate_cur': (("gs://ipbes-natcap-ecoshard-data-for-publication/MSL_Map_MERGED_Global_AVISO_NoGIA_Adjust_md5_3072845759841d0b2523d00fe9518fee.tif"), False, None, 5e3, 1),
    'slr_rcp26': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/slr_rcp26_md5_7c73cf8a1bf8851878deaeee0152dcb6.tif"), False, None, 5e3, 2),
    'slr_rcp60': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/slr_rcp60_md5_99ccaf1319d665b107a9227f2bbbd8b6.tif"), False, None, 5e3, 2),
    'slr_rcp85': ((r"gs://ipbes-natcap-ecoshard-data-for-publication/slr_rcp85_md5_3db20b7e891a71e23602826179a57e4a.tif"), False, None, 5e3, 2),
}


# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ',
    stream=sys.stdout)
LOGGER = logging.getLogger('taskgraph')
LOGGER.setLevel(logging.DEBUG)

WORKSPACE_DIR = 'global_cv_workspace'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')

if __name__ == '__main__':
    for dir_path in [WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass
