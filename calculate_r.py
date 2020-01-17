"""Given a complete global_cv vector, calculate individual Rs and Rtot."""
import argparse
import bisect
import logging
import sys

from osgeo import gdal
from osgeo import ogr
import numpy


logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

HAB_FIELDS = [
    '4_500',
    '2_2000',
    'reefs',
    'mangroves_forest',
    'saltmarsh_wetland',
    'seagrass',
    'mesoamerican_barrier_reef',
    'new_caledonian_barrier_reef',
    'great_barrier_reef',
    'keys_barrier_reef',
]


def add_cv_vector_risk(cv_risk_vector_path):
    """Use existing biophysical fields in `cv_risk_vector_path to calc total R

    Parameters:
        cv_risk_vector_path (str): path to point vector that has at least
            the following fields in it:

            * surge
            * ew
            * rei
            * slr
            * relief

            Will add the following fields:
                * Rwave
                * Rwind
                * Rsurge
                * Rrelief
                * Rslr
    Returns:
        None

    """

    cv_risk_vector = gdal.OpenEx(
        cv_risk_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    cv_risk_layer = cv_risk_vector.GetLayer()

    for base_field, risk_field in [
            ('surge', 'Rsurge'), ('ew', 'Rwave'), ('rei', 'Rwind'),
            ('slr', 'Rslr'), ('relief', 'Rrelief')]:
        cv_risk_layer.CreateField(ogr.FieldDefn(risk_field, ogr.OFTReal))
        base_array = numpy.empty(shape=(cv_risk_layer.GetFeatureCount(),))
        for index, feature in enumerate(cv_risk_layer):
            base_array[index] = feature.GetField(base_field)
        hist, bin_edges = numpy.histogram(base_array, bins=5)

        cv_risk_layer.ResetReading()
        cv_risk_layer.StartTransaction()
        for feature in cv_risk_layer:
            base_val = feature.GetField(base_field)
            risk = bisect.bisect_left(bin_edges, base_val)
            if risk < 1:
                risk = 1
            elif risk > 5:
                risk = 5
            feature.SetField(risk_field, risk)
            cv_risk_layer.SetFeature(feature)
        cv_risk_layer.CommitTransaction()
        cv_risk_vector = None
        cv_risk_layer = None
    cv_risk_layer.ResetReading()
    cv_risk_layer.CreateField(ogr.FieldDefn('Rt', ogr.OFTReal))
    for hab_field in HAB_FIELDS:
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rhab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rt_nohab_%s' % hab_field, ogr.OFTReal))
    cv_risk_layer.CreateField(
        ogr.FieldDefn('Rhab' % hab_field, ogr.OFTReal))
    cv_risk_layer.CreateField(
        ogr.FieldDefn('Rt_nohab' % hab_field, ogr.OFTReal))

    cv_risk_layer.ResetReading()
    cv_risk_layer.StartTransaction()
    for feature in cv_risk_layer:
        for hab_field in HAB_FIELDS:


        exposure_index = 1.0
        for risk_field in [
                'Rgeomorphology', 'Rhab', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                'Rrelief']:
            exposure_index *= feature.GetField(risk_field)
        exposure_index = (exposure_index)*(1./7.)
    cv_risk_layer.CommitTransaction()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate total CV risk')
    parser.add_argument('cv_risk_vector_path', help='CV Risk vector path')
    args = parser.parse_args()

    add_cv_vector_risk(args.cv_risk_vector_path)
