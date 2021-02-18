import pydicom as dicom
import numpy as np
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import scipy.ndimage as scn
from collections import defaultdict
import os
import shutil
import operator
import warnings
import math
from IPython.display import display
from dicom_contour.contour import get_ct_name_dict
import uuid
from scipy.interpolate import RegularGridInterpolator

class ArrayVolume(object):
     def __init__(self, Origin = None, Pivot= None,\
                  Vector_X= None, Vector_Y= None, Vector_Z = None,\
                  Col_Spacing= None, Row_Spacing= None, Thickness = None,\
                  Array_Volume = None, Data_Type = None, x = None, y = None, z = None, interpolating_function = None):
        self.Origin = Origin
        self.Pivot = Pivot
        self.Vector_X = Vector_X
        self.Vector_Y = Vector_Y
        self.Vector_Z = Vector_Z
        self.Col_Spacing = Col_Spacing
        self.Row_Spacing = Row_Spacing
        self.Thickness = Thickness
        self.Array_Volume= Array_Volume
        self.Data_Type = Data_Type
        self.x = x
        self.y = y
        self.z = z
        self.interpolating_function = interpolating_function
        
        
def build_image_volume(folder_path):
    image_volume = ArrayVolume()
    file_list = finder_ct_image_file(folder_path)
    init_image_ds= pydicom.dcmread(file_list[0])
    col_spacing = init_image_ds.PixelSpacing[1]
    row_spacing = init_image_ds.PixelSpacing[0]
    slice_thick = init_image_ds.SliceThickness
    start_position = dicom.dcmread(file_list[0]).ImagePositionPatient
    origin = np.array(start_position)-\
              np.array([col_spacing/2, row_spacing/2, slice_thick/2])
    if finder_rt_plan_file(folder_path)=='No Plan File':
        pivot = None
    else:
        pivot = None
        
        #pivot = dicom.dcmread(finder_rt_plan_file(folder_path))\
        #    [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
    vector_x = [1,0,0]
    vector_y = [0,1,0]
    vector_z = [0,0,1]
    data_type = 'uint'+str(init_image_ds.BitsAllocated)
    image_volume_array = []
    for frame in range(0, len(file_list)):
        image_volume_array.append\
            (dicom.dcmread(file_list[frame]).pixel_array)
    image_volume_array = np.array(image_volume_array)
    
    image_volume = ArrayVolume(origin, pivot, vector_x, vector_y, vector_z,\
                        col_spacing, row_spacing, slice_thick,
                        image_volume_array, data_type)  

    return image_volume  

def build_dose_volume(folder_path):
    dose_volume = ArrayVolume()
    ds= dicom.dcmread(finder_rt_dose_file(folder_path))
    col_spacing = ds.PixelSpacing[1]
    row_spacing = ds.PixelSpacing[0]
    slice_thick = ds.GridFrameOffsetVector[1]-ds.GridFrameOffsetVector[0]
    image_position = ds.ImagePositionPatient
    origin = np.array(image_position)-\
              np.array([col_spacing/2, row_spacing/2, slice_thick/2])
    pixel_spacing = ds.PixelSpacing
    rows = ds.Rows
    columns = ds.Columns
    x = np.arange(columns)*pixel_spacing[0] + image_position[0]
    y = np.arange(rows)*pixel_spacing[1] + image_position[1]
    z = np.array(ds.GridFrameOffsetVector) + image_position[2]
    
    pivot = None
    #pivot = dicom.dcmread(finder_rt_plan_file(folder_path))\
    #        [0x300a,0xb0][0][0x300a,0x111][0][0x300a,0x12c].value
    vector_x = [1,0,0]
    vector_y = [0,1,0]
    vector_z = [0,0,1]
    dose_volume_array = np.array(ds.pixel_array)
    data_type = 'uint'+str(ds.BitsStored)
    interpolating_function = RegularGridInterpolator((z, x, y), dose_volume_array, bounds_error=False, fill_value=0)
    
    dose_volume = ArrayVolume(origin, pivot, vector_x, vector_y, vector_z,\
                        col_spacing, row_spacing, slice_thick,
                        dose_volume_array, data_type,x,y,z, interpolating_function)  

    return dose_volume        
        
        
        
        
def finder_ct_image_file(folder_path):
    """CT Image Storage files Finder
        Locate the CT Image Storage files in the patient folder.
        Args:
            folder_path: The path to the patient folder.
        Returns:
            Return a list of the path to the CT Image Storage files.
    """
    file_list = []
    ct_image_file_list = []
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if dicom.dcmread(item).SOPClassUID.name == 'CT Image Storage':
            ct_image_file_list.append(item)
    ct_image_file_list = sorted(ct_image_file_list, 
                                key=lambda s:\
                                dicom.dcmread(s).ImagePositionPatient[2])
    return ct_image_file_list

def finder_rt_plan_file(folder_path):
    """RT Plan Storage file Finder
        Locate the RT Structure Set Storage File in the patient folder.
        Args:
            folder_path: The path to the patient folder.
        Returns:
            If there is no dose plan file in the folder, returns a string. 
            Else return the path to the RT Plan Storage file.
    """
    file_list = []
    plan_file_path = 'No Plan File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if dicom.dcmread(item).SOPClassUID.name == 'RT Plan Storage':
            plan_file_path = item

    return plan_file_path

def finder_rt_structure_file(folder_path):
    """RT Structure Set Storage file Finder
        Locate the RT Structure Set Storage File in the patient folder.
        Args:
            folder_path: The path to the patient folder.
        Returns:
            If there is no structure file in the folder, returns a string. 
            Else return the path to the RT Structure Set Storage file.
    """
    file_list = []
    structure_file_path = 'No Structure File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if dicom.dcmread(item).SOPClassUID.name == \
            'RT Structure Set Storage':
            structure_file_path = item

    return structure_file_path

def finder_rt_dose_file(folder_path):
    """RT Dose Storage file Finder
        Locate the RT Structure Set Storage File in the patient folder.
        Args:
            folder_path: The path to the patient folder.
        Returns:
            If there is no dose storage file in the folder, returns a string. 
            Else return the path to the RT Dose Storage file.
    """
    file_list = []
    dose_file_path = 'No Dose File'
    for filename in os.listdir(folder_path):
        filepath = os.path.join(folder_path, filename)
        if filepath.endswith('.dcm'):
            file_list.append(filepath)
    for item in file_list:
        if dicom.dcmread(item).SOPClassUID.name == 'RT Dose Storage':
            dose_file_path = item

    return dose_file_path        

def get_dose_file(path):
    """
    Get dose file from a given path by searching for ROIContourSequence 
    inside dicom data structure.
    More information on ROIContourSequence available here:
    http://dicom.nema.org/medical/dicom/2016c/output/chtml/part03/sect_C.8.8.6.html
    
    Inputs:
            path (str): path of the the directory that has DICOM files in it, e.g. folder of a single patient
    Return:
        dose_file (str): name of the file with the contour
    """
    # handle `/` missing
    if path[-1] != '/': path += '/'
    # get .dcm contour file
    fpaths = [path + f for f in os.listdir(path) if '.dcm' in f]
    n = 0
    dose_file = None
    for fpath in fpaths:
        f = dicom.read_file(fpath)
        if f.Modality == 'RTDOSE':
            dose_file = fpath.split('/')[-1]
            n += 1
    if n > 1: warnings.warn("There are multiple dose files, returning the last one!")
    if dose_file is None: print("No contour file found in directory")
    return dose_file

def get_roi_names(contour_data):
    """
    This function will return the names of different contour data, 
    e.g. different contours from different experts and returns the name of each.
    Inputs:
        contour_data (dicom.dataset.FileDataset): contour dataset, read by dicom.read_file
    Returns:
        roi_seq_names (list): names of the 
    """
    roi_seq_names = [roi_seq.ROIName for roi_seq in list(contour_data.StructureSetROISequence)]
    return roi_seq_names
    
def coord2pixels(contour_dataset, path, dict_name):
    """
    Given a contour dataset (a DICOM class) and path that has .dcm files of
    corresponding images. This function will return img_arr and contour_arr (2d image and contour pixels)
    Inputs
        contour_dataset: DICOM dataset class that is identified as (3006, 0016)  Contour Image Sequence
        path: string that tells the path of all DICOM images
    Return
        img_arr: 2d np.array of image with pixel intensities
        contour_arr: 2d np.array of contour with 0 and 1 labels
    """

    contour_coord = contour_dataset.ContourData

    # x, y, z coordinates of the contour in mm
    x0 = contour_coord[len(contour_coord)-3]
    y0 = contour_coord[len(contour_coord)-2]
    z0 = contour_coord[len(contour_coord)-1]
    coord = []
    for i in range(0, len(contour_coord), 3):
        x = contour_coord[i]
        y = contour_coord[i+1]
        z = contour_coord[i+2]
        l = math.sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0))
        l = math.ceil(l*2)+1
        for j in range(1, l+1):
            coord.append([(x-x0)*j/l+x0, (y-y0)*j/l+y0, (z-z0)*j/l+z0])
        x0 = x
        y0 = y
        z0 = z

    # extract the image id corresponding to given countour
    # read that dicom file (assumes filename = sopinstanceuid.dcm)
    img_ID = contour_dataset.ContourImageSequence[0].ReferencedSOPInstanceUID
    
    file_name = dict_name[img_ID]
    img = dicom.read_file(file_name)
   # img = dicom.read_file(path + img_ID + '.dcm')
    img_arr = img.pixel_array

    # physical distance between the center of each pixel
    x_spacing, y_spacing = float(img.PixelSpacing[0]), float(img.PixelSpacing[1])

    # this is the center of the upper left voxel
    origin_x, origin_y, _ = img.ImagePositionPatient

    # y, x is how it's mapped
    pixel_coords = [(np.round((y - origin_y) / y_spacing), np.round((x - origin_x) / x_spacing)) for x, y, _ in coord]

    # get contour data for the image
    rows = []
    cols = []
    for i, j in list(set(pixel_coords)):
        rows.append(i)
        cols.append(j)
    contour_arr = csc_matrix((np.ones_like(rows), (rows, cols)), dtype=np.int8, shape=(img_arr.shape[0], img_arr.shape[1])).toarray()

    return img_arr, contour_arr, img_ID


def cfile2pixels(file, path, ROIContourSeq=0):
    """
    Given a contour file and path of related images return pixel arrays for contours
    and their corresponding images.
    Inputs
        file: filename of contour
        path: path that has contour and image files
        ROIContourSeq: tells which sequence of contouring to use default 0 (RTV)
    Return
        contour_iamge_arrays: A list which have pairs of img_arr and contour_arr for a given contour file
    """
    # handle `/` missing
    if path[-1] != '/': path += '/'
    
    # get the CT id to file name dict
    dict_name = get_ct_name_dict(path)
    display(dict_name)
    
    f = dicom.read_file(path + file)
    # index 0 means that we are getting RTV information
    RTV = f.ROIContourSequence[ROIContourSeq]
    # get contour datasets in a list
    contours = [contour for contour in RTV.ContourSequence]
    img_contour_arrays = [coord2pixels(cdata, path, dict_name) for cdata in contours]  # list of img_arr, contour_arr, im_id

    # debug: there are multiple contours for the same image indepently
    # sum contour arrays and generate new img_contour_arrays
    contour_dict = defaultdict(int)
    for im_arr, cntr_arr, im_id in img_contour_arrays:
        contour_dict[im_id] += cntr_arr
    image_dict = {}
    for im_arr, cntr_arr, im_id in img_contour_arrays:
        image_dict[im_id] = im_arr
    img_contour_arrays = [(image_dict[k], contour_dict[k], k) for k in image_dict]

    return img_contour_arrays


def plot2dcontour(img_arr, contour_arr, figsize=(20, 20)):
    """
    Shows 2d MR img with contour
    Inputs
        img_arr: 2d np.array image array with pixel intensities
        contour_arr: 2d np.array contour array with pixels of 1 and 0
    """

    masked_contour_arr = np.ma.masked_where(contour_arr == 0, contour_arr)
    plt.figure(figsize=figsize)
    plt.subplot(1, 2, 1)
    plt.imshow(img_arr, cmap='gray', interpolation='none')
    plt.subplot(1, 2, 2)
    plt.imshow(img_arr, cmap='gray', interpolation='none')
    plt.imshow(masked_contour_arr, cmap='cool', interpolation='none', alpha=0.7)
    plt.show()


def slice_order(path):
    """
    Takes path of directory that has the DICOM images and returns
    a ordered list that has ordered filenames
    Inputs
        path: path that has .dcm images
    Returns
        ordered_slices: ordered tuples of filename and z-position
    """
    # handle `/` missing
    if path[-1] != '/': path += '/'
    slices = []
    for s in os.listdir(path):
        try:
            f = dicom.read_file(path + '/' + s)
            f.pixel_array  # to ensure not to read contour file
            assert f.Modality != 'RTDOSE'
            slices.append(f)
        except:
            continue

    slice_dict = {s.SOPInstanceUID: s.ImagePositionPatient[-1] for s in slices}
    ordered_slices = sorted(slice_dict.items(), key=operator.itemgetter(1))
    return ordered_slices


def get_contour_dict(contour_file, path, index):
    """
    Returns a dictionary as k: img fname, v: [corresponding img_arr, corresponding contour_arr]
    Inputs:
        contour_file: .dcm contour file name
        path: path which has contour and image files
    Returns:
        contour_dict: dictionary with 2d np.arrays
    """
    # handle `/` missing
    if path[-1] != '/': path += '/'
    # img_arr, contour_arr, img_fname
    contour_list = cfile2pixels(contour_file, path, index)

    contour_dict = {}
    for img_arr, contour_arr, img_id in contour_list:
        contour_dict[img_id] = [img_arr, contour_arr]

    return contour_dict

def get_data(path, contour_file, index):
    """
    Generate image array and contour array
    Inputs:
        path (str): path of the the directory that has DICOM files in it
        contour_file: structure file
        index (int): index of the structure
    """
    images = []
    contours = []
    # handle `/` missing
    if path[-1] != '/': path += '/'
    # get contour file
    # contour_file = get_contour_file(path)
    # get slice orders
    ordered_slices = slice_order(path)
    # get contour dict
    contour_dict = get_contour_dict(contour_file, path, index)
    
    dict_name = get_ct_name_dict(path)
    
    for k,v in ordered_slices:
        # get data from contour dict
        if k in contour_dict:
            images.append(contour_dict[k][0])
            contours.append(contour_dict[k][1])
        # get data from dicom.read_file
        else:
            file_name = dict_name[k]
            img_arr = dicom.read_file(file_name).pixel_array           
   #        img_arr = dicom.read_file(path + k + '.dcm').pixel_array
            contour_arr = np.zeros_like(img_arr)
            images.append(img_arr)
            contours.append(contour_arr)

    return np.array(images), np.array(contours)

def get_mask(path, contour_file, index, filled=True):
    """
    Generate image array and contour array
    Inputs:
        path (str): path of the the directory that has DICOM files in it
        contour_file: structure file
        index (int): index of the structure
    """
    contours = []
    # handle `/` missing
    if path[-1] != '/': path += '/'
    # get slice orders
    ordered_slices = slice_order(path)
    # get contour dict
    contour_dict = get_contour_dict(contour_file, path, index)

    for k,v in ordered_slices:
        # get data from contour dict
        if k in contour_dict:
            y = contour_dict[k][1]
            y = scn.binary_fill_holes(y) if y.max() == 1 else y
            contours.append(y)
        # get data from dicom.read_file
        else:
            img_arr = dicom.read_file(path + k + '.dcm').pixel_array
            contour_arr = np.zeros_like(img_arr)
            contours.append(contour_arr)

    return np.array(contours)


def fill_contour(contour_arr):
    # get initial pixel positions
    pixel_positions = np.array([(i, j) for i, j in zip(np.where(contour_arr)[0], np.where(contour_arr)[1])])

    # LEFT TO RIGHT SCAN
    row_pixels = defaultdict(list)
    for i, j in pixel_positions:
        row_pixels[i].append((i, j))

    for i in row_pixels:
        pixels = row_pixels[i]
        j_pos = [j for i, j in pixels]
        for j in range(min(j_pos), max(j_pos)):
            row_pixels[i].append((i, j))
    pixels = []
    for k in row_pixels:
        pix = row_pixels[k]
        pixels.append(pix)
    pixels = list(set([val for sublist in pixels for val in sublist]))

    rows, cols = zip(*pixels)
    contour_arr[rows, cols] = 1

    # TOP TO BOTTOM SCAN
    pixel_positions = pixels  # new positions added
    row_pixels = defaultdict(list)
    for i, j in pixel_positions:
        row_pixels[j].append((i, j))

    for j in row_pixels:
        pixels = row_pixels[j]
        i_pos = [i for i, j in pixels]
        for i in range(min(i_pos), max(i_pos)):
            row_pixels[j].append((i, j))
    pixels = []
    for k in row_pixels:
        pix = row_pixels[k]
        pixels.append(pix)
    pixels = list(set([val for sublist in pixels for val in sublist]))
    rows, cols = zip(*pixels)
    contour_arr[rows, cols] = 1
    return contour_arr


def create_image_mask_files(path, contour_file, index, img_format='png'):
    """
    Create image and corresponding mask files under to folders '/images' and '/masks'
    in the parent directory of path.
    
    Inputs:
        path (str): path of the the directory that has DICOM files in it, e.g. folder of a single patient
        index (int): index of the desired ROISequence
        img_format (str): image format to save by, png by default
    """
    # Extract Arrays from DICOM
    X, Y = get_data(path, contour_file, index)
    Y = np.array([scn.binary_fill_holes(y) if y.max() == 1 else y for y in Y])

    # Create images and masks folders
    new_path = '/'.join(path.split('/')[:-2])
    os.makedirs(new_path + '/images/', exist_ok=True)
    os.makedirs(new_path + '/masks/', exist_ok=True)
    for i in range(len(X)):
        plt.imsave(new_path + f'/images/image_{i}.{img_format}', X[i])
        plt.imsave(new_path + f'/masks/mask_{i}.{img_format}', Y[i])
