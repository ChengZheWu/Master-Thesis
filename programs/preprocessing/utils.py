import numpy as np
import pydicom
import cv2
import scipy.ndimage
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import math

import os

def imshow(img, grey=1):
    if grey:
        plt.imshow(np.squeeze(img), cmap="gray")
        plt.axis("off")
        plt.show()
    else:
        plt.imshow(np.squeeze(img))
        plt.show()

def check_folder(path):
    exist = os.path.exists(path)
    if not exist:
        os.makedirs(path)  # 若只用 mkdir()，如果父目錄不存在，會報錯
        print(path + "建立成功")
    else:
        print(path + "目錄已存在")
        
def get_max_coord(coords, min_size=0): # 獲得該張 dicom 的 nodule max(x, y)
    cmp = min_size
    for c in coords:
        if c > cmp:
            cmp = c
    return cmp

def get_min_coord(coords, max_size=512): # 獲得該張 dicom 的 nodule min(x, y)
    cmp = max_size
    for c in coords:
        if c < cmp:
            cmp = c
    return cmp

def cut_0(string):
    for s in range(len(string)):
        if string[s] == "0":
            result = string[s+1:]
        else:
            break
    return result

def load_scan(path):
    slices = [pydicom.read_file(s) for s in path]
    slices.sort(key = lambda x: int(x.InstanceNumber))
    
    # add slice thickness to metadata
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)

    for s in slices:
        s.SliceThickness = slice_thickness
        
    return slices

def get_pixels_hu(slices):
    #讀出像素值並且儲存成 numpy 的格式
    image = np.stack([s.pixel_array for s in slices])
    image = image.astype(np.int16)

    # Set outside-of-scan pixels to 0
    # 通常intercept是 -1024, 經過計算之後空氣大約是 0
    image[image == -2000] = 0
    
    # 轉換為Hounsfield units (HU)
    for slice_number in range(len(slices)):
        
        intercept = slices[slice_number].RescaleIntercept
        slope = slices[slice_number].RescaleSlope
        
        if slope != 1:
            image[slice_number] = slope * image[slice_number].astype(np.float64)
            image[slice_number] = image[slice_number].astype(np.int16)
            
        image[slice_number] += np.int16(intercept)
        
    img = np.array(image, dtype=np.int16)
    img[img<-1024] = -1024
    img[img>3071] = 3071
    
    return img

def normalize(image):
    MIN_BOUND = -1000
    MAX_BOUND = 400
    image = (image - MIN_BOUND) / (MAX_BOUND - MIN_BOUND)
    image[image>1] = 1
    image[image<0] = 0
    return image

def resample(image, spacing, new_spacing=[1, 1, 1]):
    spacing = np.array(spacing)
    resize_factor = spacing / new_spacing
    new_real_shape = image.shape * resize_factor
    new_shape = np.round(new_real_shape)
    real_resize_factor = new_shape / image.shape
    new_spacing = spacing / real_resize_factor
    
    image = scipy.ndimage.interpolation.zoom(image, real_resize_factor, mode='nearest')
    
    return image, new_spacing

def center_coord(img):
    coord_x = []
    coord_y = []
    for i in range(len(img)):
        for j in range(len(img)):
            if img[i][j] != 0:
                coord_x.append(i)
                coord_y.append(j)
                
    mean_x = math.floor(np.mean(coord_x))
    mean_y = math.floor(np.mean(coord_y))
    
    return mean_x, mean_y

def partition3d(nodule, size=32):
    m_x, m_y = center_coord(nodule)

    max_xbound = nodule.shape[0]
    max_ybound = nodule.shape[1]
    half_size = size//2

    if m_x - half_size < 0:
        nodule = np.pad(nodule, ((half_size-m_x, 0), (0, 0)), constant_values=0)
    elif m_x + half_size > max_xbound:
        nodule = np.pad(nodule, ((0, half_size-(max_xbound-m_x)), (0, 0)), constant_values=0)

    if m_y - half_size < 0:
        nodule = np.pad(nodule, ((0, 0), (half_size-m_y, 0)), constant_values=0)
        m_y = half_size
    elif m_y + half_size > max_ybound:
        nodule = np.pad(nodule, ((0, 0), (0, half_size-(max_ybound-m_y))), constant_values=0)

    nodule = nodule[m_x-half_size: m_x+half_size, m_y-half_size: m_y+half_size]
    
    return nodule