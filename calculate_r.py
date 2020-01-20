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
    'mangroves_forest',
    'saltmarsh_wetland',
    'seagrass',
    'reefs_all',
]

REEF_FIELDS = [
    'reefs',
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
        nan_mask = numpy.isnan(base_array)
        max_val = numpy.max(base_array[~nan_mask])
        base_array[nan_mask] = max_val
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
    cv_risk_layer.ResetReading()

    cv_risk_layer.CreateField(ogr.FieldDefn('Rt', ogr.OFTReal))
    for hab_field in HAB_FIELDS:
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rhab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rnohab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rt_nohab_%s' % hab_field, ogr.OFTReal))
    cv_risk_layer.CreateField(ogr.FieldDefn('Rhab_all', ogr.OFTReal))
    cv_risk_layer.CreateField(ogr.FieldDefn('Rt_nohab_all', ogr.OFTReal))

    # reefs are special
    cv_risk_layer.CreateField(ogr.FieldDefn('reefs_all', ogr.OFTReal))

    cv_risk_layer.ResetReading()
    cv_risk_layer.StartTransaction()
    for feature in cv_risk_layer:
        reef_risk_val = min(
            [feature.GetField(reef_field) for reef_field in REEF_FIELDS])
        feature.SetField('reefs_all', reef_risk_val)

        hab_val_map = {}
        for hab_field in HAB_FIELDS:
            hab_val = feature.GetField(hab_field)
            feature.SetField('Rhab_%s' % hab_field, hab_val)
            hab_val_map[hab_field] = hab_val

            # loop through every hab field but hab_field to calc Rhab_no
            risk_diff_list = []  # for (5-rk) vals
            for sub_hab_field in HAB_FIELDS:
                if sub_hab_field != hab_field:
                    risk_diff_list.append(5-feature.GetField(sub_hab_field))

            r_nohab = 4.8 - 0.5 * numpy.sqrt(
                (1.5 * max(risk_diff_list))**2 +
                numpy.sum([x**2 for x in risk_diff_list]) -
                max(risk_diff_list)**2)
            feature.SetField('Rnohab_%s' % hab_field, r_nohab)

        # Rhab
        # loop through every hab field but hab_field to calc Rhab_no
        risk_diff_list = []  # for (5-rk) vals
        for sub_hab_field in HAB_FIELDS:
            risk_diff_list.append(5-feature.GetField(sub_hab_field))

        r_nohab = 4.8 - 0.5 * numpy.sqrt(
            (1.5 * max(risk_diff_list))**2 +
            numpy.sum([x**2 for x in risk_diff_list]) -
            max(risk_diff_list)**2)
        feature.SetField('Rhab_all' % hab_field, r_nohab)

        # Rt
        exposure_index = 1.0
        for risk_field in [
                'Rgeomorphology', 'Rhab_all', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                'Rrelief']:
            exposure_index *= feature.GetField(risk_field)
        exposure_index = (exposure_index)**(1./7.)
        feature.SetField('Rt', exposure_index)

        # Rt_nohaball
        exposure_index = 1.0
        for risk_field in [
                'Rgeomorphology', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                'Rrelief']:
            exposure_index *= feature.GetField(risk_field)
        exposure_index = (exposure_index)**(1./6.)
        feature.SetField('Rtnohab_all', exposure_index)

        for hab_field in HAB_FIELDS:
            exposure_index = 1.0
            for risk_field in [
                    'Rgeomorphology', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                    'Rrelief', 'Rnohab_%s' % hab_field]:
                exposure_index *= feature.GetField(risk_field)
            exposure_index = (exposure_index)**(1./7.)
            feature.SetField('Rtnohab_%s' % hab_field, exposure_index)
        cv_risk_layer.SetFeature(feature)
    cv_risk_layer.CommitTransaction()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate total CV risk')
    parser.add_argument('cv_risk_vector_path', help='CV Risk vector path')
    args = parser.parse_args()

    add_cv_vector_risk(args.cv_risk_vector_path)
