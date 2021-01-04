## Scope
PyCosmic is a programme to detect cosmics in single CCD exposures from astronomical instruments, specifically  
fiber-fed integral-field spectrographs. It only excepts FITS files and returns a bad pixel mask, 
a cleaned version of the image and an error image as separate FITS file or a single one with extensions.

## Requirements
PyCosmic is running with python 3. Dependencies are listed in the requirements.txt file. 

## Installation
1. Download and unpack the PyCosmic package from github in a directory of your choice 
2. Change into the directory
3. Install PyCosmic with the normal installation command`python setup.py install`
4. Optionally: Run `pytest` to check that PyCosmic is properly working with your environment (requires pytest 
   and pytest-datafiles packages to be installed)
5. You can now try to start PyCosmic from any directory with the command `PyCosmic -h`

## Usage
IMPORTANT: The input frame is expected to be BIAS-level subtracted and coverted to ADU with the gain value. If those
steps are not yet performed. Single value bias subtraction and gain correction can be performed with PyCosmic. 
Otherwise an incorrect result is obtained. 

There are two options to use PyCosmic:
1. running from the command line with "PyCosmic"
2. loading the PyCosmic module into a customize python script

Example for both options are provided in the example director either run
* `bash run_example.sh` or 
* `python run_example.py`
and have a look at the corresponding files. 

Running `PyCosmic -h` reports any required and optimal parameter of the algorithm together with a short explanation. 

The most important parameters which control the performance of the detection algorithms are

**--siglim** The significance level of a cosmic in units of the expected noise of the pixel. 
A value of 5 should be appropriate for most applications and is the default.

**--fwhm**    The width of the Gaussian smoothing kernel along x and y dimension in pixel units. 
Prefereed is a circular case and the fwhm should be equal or smaller than the width of the instrumental PSF.

**--rlim**   Contrast threshold between the signal of cosmics and object signal in the smoothed image. 
The optimal value is not independent and coupled to the choice of --fwhm with respect to the instrumental PSF. 
For details you should consult the presentation article of PyCosmic.

**--iter**  Number of iterations to be performed by the algorithm. The default is 4 iterations and should be 
at least 3 iterations.

**--replacebox** Size of the running median box to smooth the image from potential outlier values along x and y 
direction. Should be elongated along the dispersion direction and at least 2 pixels in the cross-dispersion direction.
For other data than fiber-fed IFU the box should be symmetric with recommended 5 5 pixel size.

Additional parameter may be used as described in the programme.

## Acknowledge PyCosmic
Please cite **Husemann et al. 2012, A&A, 545, 137** if you use PyCosmic. 

### Hint on choosing fhwm and rlim
The parameters for --fwhm and --rlim need to be set properly for the given dataset to achieve good performance 
and they are dependent on each other as well as on the instrumental resolution. The table below provides some 
guidelines values, but they should be adjusted for a specific data set. 

| inst resolution | fwhm    | rlim |
| -------------- | ------- | ---- |
|  1.5     |   1.5   |  2.2 |
|  2.0     |   1.5   |  0.8	|		
|  2.0     |   2.0   |  1.0	|				
|  2.5     |   1.5   |  0.5 |
|  2.5     |   2.0   |  0.8 |
|  3.0     |   1.5   |  0.2 |
|  3.0     |   2.5   |  0.5 |

More details are found in the article.







