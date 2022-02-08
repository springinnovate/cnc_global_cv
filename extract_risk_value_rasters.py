"""Extract the risk/value rasters from the hash."""
import glob
import os
import shutil
import zipfile


if __name__ == '__main__':
    hashes = set(open('hashtrans.txt', 'r').read().rstrip().split('\n'))
    results_dir = 'results_dir'
    os.makedirs(results_dir, exist_ok=True)
    for hashline in hashes:
        root_name, short_hash = hashline.split(' -> ')
        print(root_name)
        print(short_hash)
        basename = root_name.split('_md5')[0]
        print(basename)

        value_rasters = list(glob.glob(os.path.join(
            'global_cv_workspace', short_hash, 'value_rasters', '*.tif')))
        zip_path = os.path.join(f'{basename}_value_rasters.zip')
        with zipfile.ZipFile(zip_path, 'w') as myzip:
            for value_raster in value_rasters:
                arcname = f'{basename}_{os.path.basename(value_raster)}'
                myzip.write(value_raster, arcname=arcname)

        gpkg_path = next(iter(glob.glob(
            os.path.join('global_cv_workspace', short_hash, '*.gpkg'))))
        target_gpkg_path = os.path.join(
            results_dir, f'{root_name}_point_risk.gpkg')
        shutil.copy(gpkg_path, target_gpkg_path)
