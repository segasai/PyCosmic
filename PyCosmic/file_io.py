import sys
from astropy.io import fits
import numpy as np
import argparse

__author__ = "Bernd Husemann"
__credit__ = ['Bernd Husemann', 'Sebastian Kamann', 'Christer Sandin']
__copyright__ = "Copyright 2020, Bernd Husemann"
__license__ = "MIT"
__url__ = 'https://github.com/brandherd/PyCosmic'
__maintainer__ = "Bernd Husemann"
__email__ = "berndhusemann@gmx.de"
__status__ = "Production"
__version__ = "0.6"


def read_parameters():
    parser = argparse.ArgumentParser(
        description="""
        PyCosmic programm to detect cosmic ray hits in single exposure CCD frames.  It is assumed that the input image is 
        bias subtracted and in units of electrons. Alternatively, a bias level can be subtracted and/or the image can be 
        converted from ADU to electrons with an optional gain factor. The output of the routine is a bad pixel mask and a 
        cleaned image and the internally computed Poisson+read-noise error maps.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog='PyCosmic')

    parser.add_argument("input",
                        type=str,
                        help="""File name of the input FITS file.""")
    parser.add_argument(
        "outprefix",
        type=str,
        help="""Prefix for FITS file(s) in which the output is being stored."""
    )
    parser.add_argument(
        "extension",
        type=str,
        help="""Extension number or extension name of the FITS file in which the 
        image data is stored for the procss. In case header keywords are needed for rdnoise or gain parameters,
        they need to be available in the same FITS extension.""")
    parser.add_argument(
        "rdnoise",
        type=str,
        help="""Header keyword of the CCD read-out noise in electrons 
            or alternatively the corresponding value as float number.""")
    parser.add_argument("--siglim",
                        type=float,
                        default=5.0,
                        help="""Threshold value for the significance level of 
        cosmics ray affacted pixels in units of the pixel noise.""")
    parser.add_argument("--fwhm",
                        type=float,
                        nargs=2,
                        default=[2.0, 2.0],
                        help="""FWHM in pixels of the Gaussian 
        convolution kernel used for the detection of cosmics. This should be slightly less but close to the actual 
        FWHM of the instrumental PSF in x and y direction. As a standard case we assume a symmetric PSF, but asymmetric 
        cases are also possible depending on the spectrograph and binning of the frames.  The default is 2.0,2.0."""
                        )
    parser.add_argument("--rlim",
                        type=float,
                        default=1.2,
                        help="""Threshold for the contrast value of cosmics
        against the smooth signal of the object. The optimal value depends on the FWHM of the instrument PSF and the
        FWHM of the Gaussian smoothing kernel, i.e. --fwhm. More details are provided in the publication of
        Husemann et al. 2012, AA, Volume 545, A137.""")
    parser.add_argument(
        "--iter",
        type=int,
        default=4,
        help="""Number of iteration to be performed by the algorithms. 
        Usually 4-5 iterations are needed to converge to a stable solution. Default is 5."""
    )
    parser.add_argument("--replacebox",
                        type=int,
                        nargs=2,
                        default=[5, 5],
                        help="""Size of the subimage along x and 
        y-axis around a detected cosmic ray pixel used to estimate a median value from the unaffacted pixel in order 
        to produce smooth error maps and cleaned frames. In case of fiber-fed spectrographs the box should be elongated 
        in dispersion direction and at least 2 pixels in cross-dispersion direction. For other type of data it can
        be symmetric.""")
    parser.add_argument(
        "--radius",
        type=int,
        default=0,
        help="""Integer number of neighboring pixel used to increase 
        the boundaries of the detected cosmics. Should only be used if considered to be absolutely necessary."""
    )
    parser.add_argument("--gain",
                        type=str,
                        default='1.0',
                        help="""Header keyword of the CCD gain or 
        alternatively the corresponding value as float number when the bias-subtracted image was not yet converted 
        to electron.""")
    parser.add_argument("--bias",
                        type=float,
                        default=0.0,
                        help="""Optional subtraction of a bias level. 
        The default is 0.0""")
    parser.add_argument("--ext_out",
                        action="store_true",
                        default=False,
                        help="""Flag to store the results in a single
        FITS file in extensions DATA, ERROR and MASK.""")
    parser.add_argument("--verbose",
                        action="store_true",
                        default=False,
                        help="""Flag to print some progress 
        information on the screen.""")

    args = parser.parse_args()
    input_image = args.input
    outprefix = args.outprefix
    extension = args.extension
    rdnoise = args.rdnoise
    siglim = args.siglim
    rlim = args.rlim
    fwhm = args.fwhm
    iterations = args.iter
    replace_box = args.replacebox
    radius = args.radius
    gain = args.gain
    bias = args.bias
    ext_out = args.ext_out
    verbose = args.verbose
    return input_image, outprefix, extension, rdnoise, siglim, rlim, fwhm, iterations, replace_box, radius, gain, \
           bias, ext_out, verbose


def load_file(input_image, extension_arg=0, gain_arg=1.0, rdnoise_arg=1.0):
    try:
        hdu = fits.open(input_image)
    except IOError:
        print(
            "Input FITS file could not be found. Please check file name and path."
        )
        sys.exit(1)

    try:
        try:
            extension = int(extension_arg)
        except ValueError:
            extension = extension_arg
        data = hdu[extension].data
        header_prime = hdu[0].header
        header = hdu[extension].header
    except IndexError:
        print(
            "Provided integer or string is invalid for an extension of the given FITS file."
        )
        sys.exit(2)

    try:
        gain = float(gain_arg)
    except ValueError:
        try:
            gain = float(header[gain_arg])
        except (KeyError, ValueError):
            print("Invalid header keyword provided for the gain value.")
            sys.exit(3)

    try:
        rdnoise = float(rdnoise_arg)
    except ValueError:
        try:
            rdnoise = float(header[rdnoise_arg])
        except (KeyError, ValueError):
            print("Invalid header keyword provided for the read-noise value.")
            sys.exit(4)
    return data, header, header_prime, gain, rdnoise


def save_results(result_img, outprefix, header_prime, header, flag_ext=False):

    if flag_ext:
        hdus = fits.HDUList([
            fits.PrimaryHDU,
            fits.ImageHDU(results_img.data, name='DATA'),
            fits.ImageHDU(results_img.error, name='ERROR'),
            fits.ImageHDU(result_img.mask.astype(np.uint16), name='MASK')
        ])
        hdus[0].header = header_prime
        hdus[1].header = header
        try:
            hdus.writeto(outprefix + ".pycosmic.fits",
                         output_verify='fix',
                         overwrite=True)
        except IOError:
            print(
                "Output FITS file %s could not be stored. Please check path and prefix."
                % outprefix + ".pycosmic.fits")
            sys.exit(5)
    else:
        hdu = fits.PrimaryHDU(result_img.data)
        hdu.header = header
        try:
            hdu.writeto(outprefix + ".data.fits",
                        output_verify='silentfix',
                        overwrite=True)
        except:
            print(
                "Output FITS file %s could not be stored. Please check path and prefix."
                % outprefix + ".data.fits")
            sys.exit(7)

        hdu = fits.PrimaryHDU(result_img.error)
        hdu.header = header
        hdu.writeto(outprefix + ".error.fits",
                    output_verify='silentfix',
                    overwrite=True)

        hdu = fits.PrimaryHDU(result_img.mask.astype(np.uint16))
        hdu.header = header
        hdu.writeto(outprefix + ".mask.fits",
                    output_verify='silentfix',
                    overwrite=True)
