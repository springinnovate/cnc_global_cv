"""Global CV Analysis for CNC.

Design doc is here:

https://docs.google.com/document/d/18AcJM-rXeIYgEsmqlaUwdtm7gdWLiaD6kkRkpARILlw/edit#heading=h.bbujb61ete53

"""
import logging
import os
import sys

from osgeo import gdal

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
