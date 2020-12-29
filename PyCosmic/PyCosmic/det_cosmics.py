import numpy as np
from PyCosmic.image import *

__author__ = "Bernd Husemann"
__credit__ = ['Bernd Husemann', 'Sebastian Kamann', 'Christer Sandin']
__copyright__ = "Copyright 2020, Bernd Husemann"
__license__ = "MIT"
__url__ = 'https://github.com/brandherd/PyCosmic'
__maintainer__ = "Bernd Husemann"
__email__ = "berndhusemann@gmx.de"
__status__ = "Development"
__version__ = "0.6"

def det_cosmics(data, sigma_det=5, rlim=1.2, iterations=5, fwhm_gauss=2.0, replace_box=[5, 5],
           replace_error=1e10, increase_radius=0, gain=1.0, rdnoise=1.0, bias=0.0, verbose=False):
    """
           Detects and removes cosmics from astronomical images based on Laplacian edge
           detection scheme combined with a PSF convolution approach (Husemann  et al. in prep.).

           IMPORTANT:
           The image and the readout noise are assumed to be in units of electrons.
           The image also needs to be BIAS subtracted! The gain can be entered to convert the image from ADUs to
           electros, when this is down already set gain=1.0 as the default.

            Parameters
            --------------
            image: string
                    Name of the FITS file for which the comsics should be detected
            out_mask: string
                    Name of the  FITS file with the bad pixel mask
            out_clean: string
                    Name of the  FITS file with the cleaned image
            rdnoise: float or string of header keyword
                    Value or FITS header keyword for the readout noise in electrons
            sigma_det: float, optional  with default: 5.0
                    Detection limit of edge pixel above the noise in (sigma units) to be detected as comiscs
            rlim: float, optional  with default: 1.2
                    Detection threshold between Laplacian edged and Gaussian smoothed image
            iter: integer, optional with default: 5
                    Number of iterations. Should be >1 to fully detect extended cosmics
            fwhm_gauss: float, optional with default: 2.0
                    FWHM of the Gaussian smoothing kernel in x and y direction on the CCD
            replace_box: array of two integers, optional with default: [5,5]
                    median box size in x and y to estimate replacement values from valid pixels
            replace_error: float, optional with default: 1e10
                    Error value for bad pixels in the comupted error image, will be ignored if empty
            increase_radius: integer, optional with default: 0
                    Increase the boundary of each detected cosmic ray pixel by the given number of pixels.
            verbose: bollean, optional  with default: True
                    Show information during the processing on the command line (0 - no, 1 - yes)


            References
            --------------
            B. Husemann et al. 2012  "", A&A, ??, ???

    """

    # convert all parameters to proper type
    sigma = fwhm_gauss/2.354
    box_x = int(replace_box[0])
    box_y = int(replace_box[1])

    # define Laplacian convolution kernal
    LA_kernel = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]])/4.0

    # Initiate image instances
    img_original = Image(data=data)
    img = Image(data=data)

    # subtract bias if applicable
    if (bias > 0.0) and verbose:
        print('Subtract bias level %f from image' % (bias))
    img = img - bias
    img_original = img_original - bias

    # apply gain factor to data if applicable
    if (gain != 1.0) and verbose:
        print('Convert image from ADUs to electrons using a gain factor of %f' % (gain))
    img = img * gain
    img_original = img_original * gain

    # compute noise using read-noise value
    if (rdnoise > 0.0) and verbose:
        print('A value of %f is used for the electron read-out noise.' % rdnoise)
    img_original.error = np.sqrt((np.clip(img_original.data, a_min=0.0, a_max=None) + rdnoise**2))

    select = numpy.zeros(img.dim, dtype=numpy.bool)
    img_original.mask = np.zeros(img.dim, dtype=numpy.bool)
    img.mask = np.zeros(img.dim, dtype=numpy.bool)

    # start iteration
    if verbose:
        print('Start the detection process using.')

    out = img

    for i in range(iterations):
        if verbose:
            print('Start iteration %i' % (i+1))

        # create smoothed noise fromimage
        noise = out.medianImg((box_y, box_x))
        select_neg2 = noise.data <= 0
        noise.replace_subselect(select_neg2, data=0)
        noise = (noise + rdnoise ** 2).sqrt()

        sub = img.subsample()  # subsample image
        conv = sub.convolve(LA_kernel)  # convolve subsampled image with kernel
        select_neg = conv < 0
        conv.replace_subselect(select_neg, data=0)  # replace all negative values with 0
        Lap = conv.rebin(2, 2)  # rebin the data to original resolution
        S = Lap/(noise*2)  # normalize Laplacian image by the noise
        S_prime = S-S.medianImg((5, 5))  # cleaning of the normalized Laplacian image

        # Perform additional clean using a 2D Gaussian smoothing kernel
        fine = out.convolve_gauss(sigma, sigma, mask=True)  # convolve image with a 2D Gaussian
        fine_norm = out/fine
        select_neg = fine_norm < 0
        fine_norm.replace_subselect(select_neg, data=0)
        sub_norm = fine_norm.subsample()  # subsample image
        Lap2 = sub_norm.convolve(LA_kernel)
        Lap2 = Lap2.rebin(2, 2)  # rebin the data to original resolution

        select = numpy.logical_or(numpy.logical_and(Lap2 > rlim, S_prime > sigma_det), select)

        if verbose:
            dim = img_original.dim
            det_pix = numpy.sum(select)
            print('Total number of detected cosmics: %i out of %i pixels' % (int(det_pix), dim[0] * dim[1]))

        if i == iterations-1:
            img_original.replace_subselect(select, mask=True)  # set the new mask
            if increase_radius > 0:
                mask_img = Image(data=img_original._mask)
                mask_new = mask_img.convolve(kernel=numpy.ones((2*increase_radius+1, 2*increase_radius+1)))
                img_original.mask = mask_newdata
            # replace possible corrput pixel with zeros for final output
            out = img_original.replaceMaskMedian(box_x, box_y, replace_error=replace_error)
        else:
            out.replace_subselect(select, mask=True)  # set the new mask
            out = out.replaceMaskMedian(box_x, box_y, replace_error=None)  # replace possible corrput pixel with zeros

    return out