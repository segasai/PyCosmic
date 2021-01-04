from astropy.io import fits as pyfits
import numpy
from scipy import ndimage
from scipy import stats

__author__ = "Bernd Husemann"
__credit__ = ['Bernd Husemann', 'Sebastian Kamann', 'Christer Sandin']
__copyright__ = "Copyright 2020, Bernd Husemann"
__license__ = "MIT"
__url__ = 'https://github.com/brandherd/PyCosmic'
__maintainer__ = "Bernd Husemann"
__email__ = "berndhusemann@gmx.de"
__status__ = "Production"
__version__ = "0.6"

class Image(object):
    def __init__(self, data=None, error=None, mask=None):
        self._data = data
        if self._data is not None:
            self._dim = self._data.shape
        else:
            self._dim = None
        self._mask = mask
        self._error = error
    
    def __add__(self, other):
        """
        Operator to add two Images or add another type if possible
        """
        if isinstance(other, Image):
            # define behaviour if the other is of the same instance
   
            img = Image()
            
            # add data if contained in both
            if self._data is not None and other._data is not None:
                new_data = self._data+other._data 
                img.data = new_data
            else:
                img.data = self._data
            
            # add error if contained in both 
            if self._error is not None and other._error is not None:
                new_error = numpy.sqrt(self._error**2+other._error**2) 
                img.error = new_error
            else:
                img.error = self._error
                
            # combined mask of valid pixels if contained in both     
            if self._mask is not None and other._mask is not None:
                new_mask = numpy.logical_or(self._mask, other._mask) 
                img.mask = new_mask
            else:
                img.mask = self._mask
            return img

        elif isinstance(other,  numpy.ndarray):
            img = Image(error=self._error, mask=self._mask)
    
            if self._data is not None:  # check if there is data in the object
                dim = other.shape
                # add ndarray according do its dimensions
                if self._dim == dim:
                    new_data = self._data+other
                elif len(dim) == 1:
                    if self._dim[0] == dim[0]:
                        new_data = self._data+other[:, numpy.newaxis]
                    elif self._dim[1] == dim[0]:
                        new_data = self._data+other[numpy.newaxis, :]
                else:
                    new_data = self._data
                img.data = new_data
            return img
        else:
            # try to do addition for other types, e.g. float, int, etc.
            try:
                new_data = self._data+other
                img = Image(data=new_data, error=self._error, mask=self._mask)
                return img
            except:
                # raise exception if the type are not matching in general
                raise TypeError("unsupported operand type(s) for +: %s and %s" %
                                (str(type(self)).split("'")[1], str(type(other)).split("'")[1]))
                
    def __radd__(self, other):
        self.__add__(other)
        
    def __sub__(self, other):
        """
        Operator to subtract two Images or subtract another type if possible
        """
        if isinstance(other, Image):
            # define behaviour if the other is of the same instance
   
            img = Image()
            
            # subtract data if contained in both
            if self._data is not None and other._data is not None:
                new_data = self._data-other._data 
                img.data = new_data
            else:
                img.data = self._data
            
            # add error if contained in both 
            if self._error is not None and other._error is not None:
                new_error = numpy.sqrt(self._error**2+other._error**2) 
                img.error = new_error
            else:
                img.error = self._error
                
            # combined mask of valid pixels if contained in both     
            if self._mask is not None and other._mask is not None:
                new_mask = numpy.logical_or(self._mask, other._mask) 
                img.mask = new_mask
            else:
                img.mask = self._mask
            return img

        elif isinstance(other,  numpy.ndarray):
            img = Image(error=self._error, mask=self._mask)
    
            if self._data is not None:  # check if there is data in the object
                dim = other.shape
                # add ndarray according do its dimensions
                if self._dim == dim:
                    new_data = self._data-other
                elif len(dim) == 1:
                    if self._dim[0] == dim[0]:
                        new_data = self._data-other[:, numpy.newaxis]
                    elif self._dim[1] == dim[0]:
                        new_data = self._data-other[numpy.newaxis, :]
                else:
                    new_data = self._data
                img.data = new_data
            return img
        else:
            # try to do addtion for other types, e.g. float, int, etc.
            try:
                new_data = self._data-other
                img = Image(data=new_data, error=self._error, mask=self._mask)
                return img
            except:
                # raise exception if the type are not matching in general
                raise TypeError("unsupported operand type(s) for -: %s and %s" %
                                (str(type(self)).split("'")[1], str(type(other)).split("'")[1]))

    def __truediv__(self, other):
        """
        Operator to divide two Images or divide by another type if possible
        """
        if isinstance(other, Image):
            # define behaviour if the other is of the same instance
   
            img = Image()
            
            # subtract data if contained in both
            if self._data is not None and other._data is not None:
                new_data = self._data/other._data 
                img.data = new_data
            else:
                img.data = self._data
            
            # add error if contained in both 
            if self._error is not None and other._error is not None:
                new_error = numpy.sqrt((self._error/other._data)**2+(self._data*other._error/other._data**2)**2) 
                img.error = new_error
            else:
                img.error = self._error
                
            # combined mask of valid pixels if contained in both     
            if self._mask is not None and other._mask is not None:
                new_mask = numpy.logical_or(self._mask, other._mask) 
                img.mask = new_mask
            else:
                img.mask = self._mask
            return img

        elif isinstance(other,  numpy.ndarray):
            img = Image(error=self._error, mask=self._mask)
    
            if self._data is not None:  # check if there is data in the object
                dim = other.shape
                # add ndarray according do its dimensions
                if self._dim == dim:
                    new_data = self._data/other
                    if self._error is not None:
                        new_error = self._error/other
                    else:
                        new_error = None
                elif len(dim) == 1:
                    if self._dim[0] == dim[0]:
                        new_data = self._data/other[:, numpy.newaxis]
                        if self._error is not None:
                            new_error = self._error/other[:, numpy.newaxis]
                        else:
                            new_error = None
                    elif self._dim[1] == dim[0]:
                        new_data = self._data/other[numpy.newaxis, :]
                        if self._error is not None:
                            new_error = self._error/other[numpy.newaxis, :]
                        else:
                            new_error = None
                else:
                    new_data = self._data
                img.data = new_data
                img.error = new_error
            return img
        else:
            # try to do addtion for other types, e.g. float, int, etc.
            try:
                new_data = self._data/other
                if self._error is not None:
                    new_error = self._error/other
                else:
                    new_error = None
                img = Image(data=new_data, error=new_error, mask=self._mask)
                return img
            except:
                # raise exception if the type are not matching in general
                raise TypeError("unsupported operand type(s) for /: %s and %s" %
                                (str(type(self)).split("'")[1], str(type(other)).split("'")[1]))
                
    def __mul__(self, other):
        """
        Operator to divide two Images or divide by another type if possible
        """
        if isinstance(other, Image):
            # define behaviour if the other is of the same instance
   
            img = Image()
            
            # subtract data if contained in both
            if self._data is not None and other._data is not None:
                new_data = self._data*other._data 
                img.data = new_data
            else:
                img.data = self._data
            
            # add error if contained in both 
            if self._error is not None and other._error is not None:
                new_error = numpy.sqrt((self._error*other._data)**2+(self._data*other._error)**2) 
                img.error = new_error
            else:
                img.error = self._error
                
            # combined mask of valid pixels if contained in both     
            if self._mask is not None and other._mask is not None:
                new_mask = numpy.logical_or(self._mask, other._mask) 
                img.mask = new_mask
            else:
                img.mask = self._mask
            return img

        elif isinstance(other,  numpy.ndarray):
            img = Image(error=self._error, mask=self._mask)
    
            if self._data is not None:  # check if there is data in the object
                dim = other.shape
                # add ndarray according do its dimensions
                if self._dim == dim:
                    new_data = self._data*other
                elif len(dim) == 1:
                    if self._dim[0] == dim[0]:
                        new_data = self._data*other[:, numpy.newaxis]
                    elif self._dim[1] == dim[0]:
                        new_data = self._data*other[numpy.newaxis, :]
                else:
                    new_data = self._data
                img.data = new_data
            return img
        else:
            # try to do addtion for other types, e.g. float, int, etc.
            try:
                new_data = self._data*other
                img = Image(data=new_data, error=self._error, mask=self._mask)
                return img
            except:
                # raise exception if the type are not matching in general
                raise TypeError("unsupported operand type(s) for *: %s and %s" %
                                (str(type(self)).split("'")[1], str(type(other)).split("'")[1]))

    ## define comparison operators as a comparison with
    def __rmul__(self, other):
         self.__mul__(other)

    def __lt__(self, other):
         return self._data < other

    def __le__(self, other):
         return self._data <= other

    def __eq__(self, other):
         return self._data == other

    def __ne__(self, other):
         return self._data != other

    def __gt__(self, other):
         return self._data > other

    def __ge__(self, other):
         return self._data >= other

    def sqrt(self):
        """
            Computes the square root  of the image 
            
            Returns
            -----------
            Image : data_model.Image object
                A full Image object
                       
        """
        if self._data is not None:
            new_data = numpy.sqrt(self._data)   # sqrt of the data
        else:
            new_data = None
        
        if self._error is not None and self._data is not None:
            new_error = 1/(2*new_data)*self._error   # corresponding error
        else:
            new_error = None
        #   return new Image object with corresponding data
        return Image(data=new_data, error=new_error,  mask=self._mask)

    @property
    def dim(self):
        """
            Returns the dimension of the image 
            
            Returns
            -----------
            _dim :  tuple 
                The dimension of the image (y,x)
                       
        """
        return self._dim

    @property
    def data(self):
        """
            Returns the stored data of the image

            Returns
            -----------
            _data :  numpy.ndarray
            The stored data of the image

        """
        return self._data

    @data.setter
    def data(self, data):
        if data is None:
            self._data = None
        else:
            if self._dim is None:
                self._data = data
                self._dim = data.shape
            elif data.shape == self._dim:
                self._data = data
            else:
                raise RuntimeError("Incompatible dimensions for replacing Image data.")

    @property
    def mask(self):
        """
            Returns the bad pixel mask of the image 
            
            Returns
            -----------
            _mask :  numpy.ndarray
                The bad pixel mask of the image
                       
        """
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = None
        else:
            if self._dim is None:
                self._mask = mask
                self._dim = mask.shape
            elif mask.shape == self._dim:
                self._mask = mask
            else:
                raise RuntimeError("Incompatible dimensions for replacing Image mask.")

    @property
    def error(self):
        """
            Returns the associated error of the image 
            
            Returns
            -----------
            _error :  numpy.ndarray
                The associated error of the image
                       
        """
        return self._error

    @error.setter
    def error(self, error):
        if error is None:
            self._error = None
        else:
            if self._dim is None:
                self._error = error
                self._dim = error.shape
            elif error.shape == self._dim:
                self._error = error
            else:
                raise RuntimeError("Incompatible dimensions for replacing Image error.")

    def replace_subselect(self, select, data=None, error=None, mask=None):
        """
            Set data for an Image. Specific data values can replaced according to a specific selection.

            Parameters
            --------------
            select : numpy.ndarray(bool)
                array defining the selection of pixel to be set
            data : numpy.ndarray(float), optional with default = None
                array corresponding to the data to be set
            error : numpy.ndarray(float), optional with default = None
                array corresponding to the data to be set
            mask : numpy.ndarray(bool), optional with default = None
                array corresponding to the bad pixel to be set
        """
        if data is not None:
            self._data[select] = data
        if mask is not None:
            self._mask[select] = mask
        if error is not None:
            self._error[select] = error

    def replaceMaskMedian(self, box_x, box_y, replace_error=1e20):
        """
            Replace bad pixels with the median value of pixel in a rectangular filter window 
            
            Parameters
            --------------
            box_x : int
                Pixel size of filter window in x direction
            box_y : int 
                Pixel size of filter window in y direction
            replace_error : float, optional with default: None
                Error that should be set for bad pixel
                
            Returns
            -----------
            new_image :  Image object
                Subsampled image
        """

        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")

        idx = numpy.indices(self._dim)  # create an index array
        # get x and y coordinates of bad pixels
        
        y_cors = idx[0][self._mask]
        x_cors = idx[1][self._mask]
        
        out_data = self._data
        out_error = self._error
        
        # esimate the pixel distance form the bad pixel to the filter window boundary
        delta_x = numpy.ceil(box_x/2.0)
        delta_y = numpy.ceil(box_y/2.0)
        
        # iterate over bad pixels
        for m in range(len(y_cors)):
            # computes the min and max pixels of the filter window in x and y 
            range_y = numpy.clip([y_cors[m]-delta_y, y_cors[m]+delta_y+1], 0, self._dim[0]-1).astype(numpy.uint16)
            range_x = (numpy.clip([x_cors[m]-delta_x, x_cors[m]+delta_x+1], 0, self._dim[1]-1)).astype(numpy.uint16)
            # compute the masked median within the filter window and replace data
            select = self._mask[range_y[0]:range_y[1], range_x[0]:range_x[1]] == 0
            out_data[y_cors[m], x_cors[m]] = numpy.median(self._data[range_y[0]:range_y[1],
                                                          range_x[0]:range_x[1]][select])
            if self._error is not None and replace_error is not None:
                # replace the error of bad pixel if defined
                out_error[y_cors[m], x_cors[m]] = replace_error
                
        # create new Image object
        new_image = Image(data=out_data, error=out_error,  mask=self._mask)
        return new_image

    def subsample(self):
        """
            Subsample the image by a factor of 2, e.g. each pixel is divided into 4 pixel so that their sum is 4
            times the original one.
            
            Returns
            -----------
            new_image :  Image object
                Subsampled image
                       
        """
        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")

        # create empty array with 2 time larger size in both axes
        new_dim = (self._dim[0]*2, self._dim[1]*2)
        new_data = numpy.zeros(new_dim, dtype=numpy.float32)
        if self._error is not None:
            new_error = numpy.zeros(new_dim, dtype=numpy.float32)
        else:
            new_error = None
        if self._mask is not None:
            new_mask = numpy.zeros(new_dim, dtype='bool')
        else:
            new_mask = None

        # set pixel for the subsampled data, error and mask
        new_data[::2, ::2] = self._data
        new_data[::2, 1::2] = self._data
        new_data[1::2, ::2] = self._data
        new_data[1::2, 1::2] = self._data
        if self._error is not None:
            new_error[::2, ::2] = self._error
            new_error[::2, 1::2] = self._error
            new_error[1::2, ::2] = self._error
            new_error[1::2, 1::2] = self._error
        if self._mask is not None:
            new_mask[::2, ::2] = self._mask
            new_mask[::2, 1::2] = self._mask
            new_mask[1::2, ::2] = self._mask
            new_mask[1::2, 1::2] = self._mask
        
        # create new Image object with the new subsample data    
        new_image = Image(data=new_data, error=new_error,  mask=new_mask)
        return new_image      

    def rebin(self, bin_x, bin_y):
        """
            Rebin the image by regullarly summing up the pixel in a regual rectangular binning window with size bin_x
            times bin_y. Make sure that the size of the binning window matches with the total number of pixel in
            the original image.
            
            Parameters
            --------------
            bin_x : int
                Pixel size of the binning window in x direction
            bin_y : int 
                Pixel size of the binning window in y direction
                
            Returns
            -----------
            new_image :  Image object
                Subsampled image
        """
        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")
        if (self._dim[0] % bin_y) != 0 or (self._dim[1] % bin_x) != 0:
            raise RuntimeError("Binning cannot be performed. Input dimensions are not a multiple of the requested"
                               " binning.")

        # sum over the data array over each axis by the given pixel 
        new = numpy.sum(numpy.reshape(self._data, (self._dim[0], int(self._dim[1]/bin_x), int(bin_x))), 2)
        new2 = numpy.sum(numpy.reshape(new, (int(self._dim[0]/bin_y), int(bin_y), int(self._dim[1]/bin_x))), 1)
        
        if self._error is not None:
            # sum over the error array (converted to variance and back) over each axis by the given pixel
            error_new = numpy.sum(numpy.reshape(self._error**2, (self._dim[0], int(self._dim[1]/bin_x), int(bin_x))), 2)
            error_new2 = numpy.sqrt(numpy.sum(numpy.reshape(error_new, (int(self._dim[0]/bin_y), int(bin_y),
                                                                        int(self._dim[1]/bin_x))), 1))
        else:
            error_new2 = None
            
        if self._mask is not None:
            # create the new  bad pixel mask 
            mask_new = numpy.sum(numpy.reshape(self._mask, (self._dim[0], int(self._dim[1]/bin_x), int(bin_x))), 2)
            mask_new2 = numpy.sum(numpy.reshape(mask_new, (int(self._dim[0]/bin_y), int(bin_y), int(self._dim[1]/bin_x))), 1)
            # if only one bad pixel in the binning pixel exists the binned pixel will have the bad pixel status
            new_mask = mask_new2 > 0
        else:
            new_mask = None
        # create new Image object and return
        new_img = Image(data=new2, error=error_new2, mask=new_mask)
        return new_img
        
    def convolve(self, kernel, mode='nearest'):
        """
            Convolves the data of the Image with a given kernel. The mask and error information will be unchanged.
            
            Parameters
            --------------
            kernel : ndarray
                Convolution kernel
            mode :  string, optional with default: 'nearest'
                Set the mode how to handle the boundarys within the convolution
  
                
            Returns
            -----------
            new_image :  Image object
                Convolved image
        """
        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")

        # convolve the data array with the given convolution kernel
        new = ndimage.filters.convolve(self._data, kernel, mode=mode)
        if self._error is not None:
            new_error = numpy.sqrt(ndimage.filters.convolve(self._error**2, kernel, mode=mode))
        else:
            new_error = None
        # create new Image object with the error and the mask unchanged and return
        new_image = Image(data=new, error=new_error,  mask=self._mask)
        return new_image
        
    def convolve_gauss(self, sigma_x, sigma_y, mode='nearest', mask=False):
        """
            Convolves the data of the Image with a given kernel. The mask and error information will be unchanged.
            
            Parameters
            --------------
            sigma_x : float
                With of the Gaussian in pixels along the x direction
            sigma_y : float
                With of the Gaussian in pixels along the y direction
            mode :  string, optional with default: 'nearest'
                Set the mode how to handle the boundarys within the convolution

                
            Returns
            -----------
            new_image :  Image object
                Convolved Image
        """

        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")

        # convolve the data array with the 2D Gaussian convolution kernel

        if self._mask is not None and mask:
            mask_data = self._data[self._mask]
            self._data[self._mask] = 0
            gauss = ndimage.filters.gaussian_filter(self._data, (sigma_y, sigma_x), mode=mode)
            scale = ndimage.filters.gaussian_filter((self._mask == False).astype('float32'), (sigma_y, sigma_x),
                                                    mode=mode)
            new = gauss/scale
            self._data[self._mask] = mask_data
        else:
            new = ndimage.filters.gaussian_filter(self._data, (sigma_y, sigma_x), mode=mode)
        # create new Image object with the error and the mask unchanged and return
        new_image = Image(data=new, error=self._error,  mask=self._mask)
        return new_image
        
    def medianImg(self, size, mode='nearest'):
        """
            Return a new Image that has been median filtered with a filter window of given size.
            
            Parameters
            --------------
            size : tuple of int
                Size of the filter window 
            mode : string, optional with default: nearest
                Set the mode how to handle the boundarys within the convolution
                Possilbe modes are: reflect, constant, nearest, mirror,  wrap
                
            Returns
            -----------
            image :  Image object
                An Image object with the median filter data
        """

        if self._data is None:
            raise RuntimeError("Image object is empty. Nothing to process.")

        # applying the median filter
        new_data = ndimage.filters.median_filter(self._data, size, mode=mode)
        # create a new Image object
        image = Image(data=new_data, error=self._error,  mask=self._mask)
        return image
