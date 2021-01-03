from astropy.io import fits
import PyCosmic
import numpy as np

### explicit loading of image header and header keywords
hdu = fits.open('PMAS_exp2.fits')
img = hdu[1].data
hdr = hdu[0].header
gain = hdr['GAIN']
rdnoise=hdr['RDNOISE']
header_prime = hdr
header = hdu[1].header
hdu.close()

### alternatively one can use a convenience function from PyCosmic
#(img, header, header_prime, gain, rdnoise) = PyCosmic.load_file('PMAS_exp2.fits', 1, gain_arg=1.2, rdnoise_arg=2.3)

#### Run the PyCosmic detection algorithm
out = PyCosmic.det_cosmics(img,gain=gain,rdnoise=rdnoise, rlim=1.4,iterations=4, replace_box=[10,2],
                  replace_error=100, verbose=True)

### explicit writing of the resulting images
hdu = fits.PrimaryHDU(out.data)
hdu.writeto('PMAS_exp2.data.fits',overwrite=True)

hdu = fits.PrimaryHDU(out.error)
hdu.writeto('PMAS_exp2.error.fits',overwrite=True)

hdu = fits.PrimaryHDU(out.mask.astype(np.uint16))
hdu.writeto('PMAS_exp2.mask.fits',overwrite=True)

### alternatively one can use a convenience function from PyCosmic
#PyCosmic.save_results(out, 'PMAS_exp2', header_prime, header, flag_ext=False)



