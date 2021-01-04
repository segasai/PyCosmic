import pytest
import os
import numpy as np
from astropy.io import fits
from scipy import ndimage
import PyCosmic


#@pytest.mark.datafiles('../example/PMAS_exp2.fits')
#def test_load(datafiles):
#     path = str(datafiles)
#     file_path = os.path.join(path, 'PMAS_exp2.fits')
#     assert load_file(input_image)


@pytest.mark.datafiles('example/PMAS_exp2.fits')
def test_detection(datafiles):
     path = str(datafiles)
     file_path = os.path.join(path, 'PMAS_exp2.fits')
     hdu = fits.open(file_path)
     img = hdu[1].data
     hdr = hdu[0].header
     out = PyCosmic.det_cosmics(img, gain=hdr['GAIN'], rdnoise=hdr['RDNOISE'], rlim=1.4, replace_box=[10, 2],
                                replace_error=100, iterations=4, verbose=False)
     assert np.sum(out.mask) == 1514
     assert np.mean(out.error[out.mask]) == 100

