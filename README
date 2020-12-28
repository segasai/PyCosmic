----------INSTALLATION--------------

A) For a Linux/Unix system with root permissions:

1. Unpack the compressed tar file in a directory of your choice (DIR)
2. run "sudo python setup.py install"
3. You can now try to start PyCosmic from any directory with the command,
      "PyCosmic"

B) For a Linux/Unix system without root permissions:

1. Unpack the compressed tar file in a directory of your choice (DIR)
2. run "sudo python setup.py install --home=HOME" 
HOME is a directory of your choice in which you have permission to write. The python module when then be install in the subdirectories "HOME/lib/python" and the executable in the directory "HOME/bin". Please consult the page http://docs.python.org/install/index.html#alternate-installation-the-home-scheme for further information.
3. make sure that the module and executable pathes are included in the PYTHONPATH shell variable and the PATH shell variable, respectively.
4. You can now try to start PyCosmic from any directory with the command,
      "PyCosmic"


----------USAGE--------------

PyCosmic is a programme to detect cosmics in single exposure of astronomical instruments, specifically  fiber-fed integral-field spectrographs. It only excepts FITS files and returns a bad pixel mask FITS image and an image cleaned from any detected cosmics. 

IMPORTANT: The input frame is expected to be BIAS subtracted in any case! Without proper BIAS subtraction the algorithm will return incorrect results.

"PyCosmic -h" will report any required and optimal parameter of the algorithm together with a short explanation. 

The most important parameters which control the performance of the detection algorithms are

--siglim    The significance level of a cosmic in units of the expected noise of the pixel. A value of 5 should be appropriate for most applications and is the default.
--fwhm    The width of the Gaussian smoothing kernel in pixel units. Should be equal or smaller than the width of the instrumental PSF.
--rlim   Contrast threshold between the signal of cosmics and object signal in the smoothed image. The optimal value is not independent and coupled to the choice of --fwhm with respect to the instrumental PSF. For details you should consult the presentation article of PyCosmic. 
--iter  Number of iterations to be performed by the algorithm. The default is 5 iterations and should be at least 3 iterations.

Additional parameter may be used as described in the programme.


In case you have used PyCosmic for reducing your data. Please cite the following reference in your article(s):
Husemann et al. 2012, A&A, 545, 137


----------EXAMPLE--------------

PyCosmic IMAGE_IN.fits MASK_OUT.fits CLEAN_OUT.fits RDNOISE --fwhm 2 --rlim 1 --iter 6

Here IMAGE_IN.fits, MASK_OUT.fits, and CLEAN_OUT.fits need to replaced by the corresponding file names. RDNOISE is either a header keyword or a float representative for the read-out noise of the CCD.

The parameters for --fwhm and --rlim need to be set properly for the given dataset to achieve good performance and they are dependent on each other as well as on the instrumental resolution. The table below provides some guidelines values 

inst. resolution , fwhm  , rlim 
    1.5 	    1.5	    2.2
    2.0		    1.5	    0.8			
    2.0		    2.0	    1.0					
    2.5		    1.5     0.5
    2.5		    2.0     0.8 
    3.0		    1.5     0.2
    3.0		    2.5     0.5

You will certainly need to adjust the value for rlim to achieve the optimal results, but you should be close already. For further details please have a look at the presentation article of PyCosmic.


----------GETTING HELP--------------
In case of troubles with the software and other questions, please consult the project page at pycosmic.sf.net to open a ticket or a discussion in the Forum. 







