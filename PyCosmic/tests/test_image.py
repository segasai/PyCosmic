import os
import pytest
import numpy as np
import PyCosmic

@pytest.fixture()
def array1():
    test_array1 = np.ones((6,6),dtype = np.float32) * 4
    test_array1[0:3,0:3] = 16.0
    return test_array1

@pytest.fixture()
def array2():
    test_array2 = np.ones((6,6),dtype = np.float32)*2
    return test_array2

@pytest.fixture()
def array3():
    test_array3 = np.ones((6,6),dtype = np.float32)
    test_array3[3,3] = -1
    return test_array3

@pytest.fixture()
def array4():
    test_array4 = np.zeros((6,6),dtype = bool)
    test_array4[3,3] = True
    return test_array4

def test_image_empty():
    img_empty = PyCosmic.Image()
    assert img_empty.data is None
    assert img_empty.dim is None
    assert img_empty.mask is None
    assert img_empty.error is None

def test_image_addition():
    img_empty = PyCosmic.Image()
    assert (img_empty + img_empty).data is None
    assert (img_empty + img_empty).dim is None
    assert (img_empty + img_empty).mask is None
    assert (img_empty + img_empty).error is None

def test_image_multiplication():
    img_empty = PyCosmic.Image()
    assert (img_empty * img_empty).data is None
    assert (img_empty * img_empty).dim is None
    assert (img_empty * img_empty).mask is None
    assert (img_empty * img_empty).error is None

def test_image_subtraction():
    img_empty = PyCosmic.Image()
    assert (img_empty - img_empty).data is None
    assert (img_empty - img_empty).dim is None
    assert (img_empty - img_empty).mask is None
    assert (img_empty - img_empty).error is None

def test_image_division():
    img_empty = PyCosmic.Image()
    assert (img_empty / img_empty).data is None
    assert (img_empty / img_empty).dim is None
    assert (img_empty / img_empty).mask is None
    assert (img_empty / img_empty).error is None

def test_image_sqrt():
    img_empty = PyCosmic.Image()
    assert img_empty.sqrt().data is None
    assert img_empty.sqrt().dim is None
    assert img_empty.sqrt().error is None


def test_image_simple(array1,array2,array3,array4):
    img1 = PyCosmic.Image(data=array1, error=array2, mask=array4)
    assert np.sum(img1.data)==252
    assert len(img1.dim) == 2
    assert np.sum(img1.dim) == 12
    assert np.sum(img1.error) == 72
    assert np.sum(img1.mask) == 1

    img_sqrt = img1.sqrt()
    assert np.sum(img_sqrt.data) == 90
    assert np.sum(img_sqrt.error) == np.sum(1.0/(2*array1**0.5)*array2)
    assert np.array_equal(img_sqrt.mask,array4)

    img_sum = img1 + img1
    assert np.sum(img_sum.data) == 252*2
    assert np.sum(img_sum.dim) == 12
    assert np.sum(img_sum.error) == np.sum(np.sqrt(2*array2**2))
    assert np.sum(img_sum.mask) == 1

def test_image_subsample(array1, array2, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        sub = img_empty.subsample()

    img1 = PyCosmic.Image(data=array1)
    sub = img1.subsample()
    assert np.sum(sub.data) == 252*4
    assert np.sum(sub.dim) == 24
    assert sub.mask is None
    assert sub.error is None

    img2 = PyCosmic.Image(data=array1, error=array2, mask=array4)
    sub = img2.subsample()
    assert np.sum(sub.data) == 252 * 4
    assert np.sum(sub.dim) == 12 * 2
    assert np.sum(sub.mask) == 1 * 4
    assert np.sum(sub.mask.shape) == 24
    assert np.sum(sub.error) == 72 * 4
    assert np.sum(sub.error.shape) == 24

def test_image_convolve(array1, array2, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        img_empty.convolve(0)

    LA_kernel = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]])/4.0
    img1 = PyCosmic.Image(data=array1)
    conv = img1.convolve(LA_kernel)
    assert np.sum(conv.data) == 0
    assert np.sum(conv.dim) == 12

    expand_kernel = np.ones((2+1, 2+1))
    conv = img1.convolve(expand_kernel)
    assert np.sum(conv.data) == 2268

    img2 = PyCosmic.Image(data=array4)
    conv = img2.convolve(expand_kernel)
    assert np.sum(conv.data) == 9

    img3 = PyCosmic.Image(data=array1, error=array2, mask=array4)
    conv = img3.convolve(LA_kernel)
    assert np.sum(conv.data) == 0
    assert np.sum(conv.dim) == 12
    assert np.sum(conv.error) == 0
    assert conv.mask is not None

def test_image_subselect(array1, array2, array4):
    select = array4
    img = PyCosmic.Image(data=array1, error=array2, mask=array4)
    img.replace_subselect(select,data=0)
    img.replace_subselect(select,error=0)
    img.replace_subselect(select, mask=False)
    assert np.sum(img.error) == 72 - 2
    assert np.sum(img.data) == 252 - 4
    assert np.sum(img.mask) == 0

def test_image_rebin(array1, array2, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        img_empty.rebin(1, 1)
    img = PyCosmic.Image(data=array1, error=array2, mask=array4)
    binned = img.rebin(2, 2)
    assert np.sum(binned.data) == 252
    assert np.sum(binned.dim) == 6
    assert np.sum(binned.error) == 4*9
    assert np.sum(binned.mask)  == 1
    with pytest.raises(RuntimeError, match="Binning cannot be performed. Input dimensions are not a multiple of the "
                                           "requested binning."):
        binned = img.rebin(5, 5)

def test_image_medianImg(array1, array2, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        img_empty.medianImg(0)
    img = PyCosmic.Image(data=array1, error=array2, mask=array4)
    med = img.medianImg([3,3],mode='nearest')
    assert np.sum(med.data) == 240.0
    assert np.sum(med.dim) == 12
    assert np.sum(med.mask) == 1
    assert np.sum(med.error) == 72
    med = img.medianImg([3, 1], mode='nearest')
    assert np.sum(med.data) == 252.0
    med = img.medianImg(2, mode='nearest')
    assert np.sum(med.data) == 324.0

def test_image_convolve_gauss(array1, array2, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        img_empty.convolve_gauss(0,0)
    img = PyCosmic.Image(data=array1, error=array2, mask=array4)
    smooth1 = img.convolve_gauss(1, 1, mode='nearest')
    assert np.sum(smooth1.data) == 252.0
    assert np.sum(smooth1.dim) == 12
    smooth2 = img.convolve_gauss(1, 1, mode='nearest', mask=True)
    assert np.sum(smooth1.data) < np.sum(smooth2.data)

def test_image_replaceMask(array3, array4):
    img_empty = PyCosmic.Image()
    with pytest.raises(RuntimeError, match="Image object is empty. Nothing to process."):
        img_empty.replaceMaskMedian(10,10)

    img = PyCosmic.Image(data=array3, mask=array4)
    assert img.data[3,3] == -1
    img.replaceMaskMedian(3,3)
    assert img.data[3, 3] == 1

