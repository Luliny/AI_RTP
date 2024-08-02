import os
import sys
import SimpleITK as sitk
from skimage import transform
from pydicom import dcmread
import numpy as np
import cv2
import argparse
import torch
#from data import create_dataset
import pdb
import re
import pandas as pd

sys.path.append('.')

from DICOM_intake import dicom_validation

def _get_tag(doc=None, group_num=None, element_num=None, default=None, string=True):
    """
    Helper function to extract specific DICOM tag value.

    Args:
        doc (pydicom object): DICOM object
        group_num: 4-digit hexadecimal number (e.g. 0x1010) corresponding to first half of DICOM tag e.g. (1010,xxxx)
        element_num: 4-digit hexadecimal number (e.g. 0x3022) corresponding to second half of DICOM tag e.g. (xxxx,3022)
        default (object): default value to return in place of any errors or missing values
        string (boolean): specifies whether to extract value as a string (default: True)
    
    Returns: 
        object: Requested content from the DICOM object
    """
    if not doc or not hasattr(doc, 'keys'):
        return(default)
    if not group_num or not element_num:
        return(default)
    if [group_num,element_num] in doc:
        doc_field = doc[group_num,element_num]
        if doc_field.VR == 'SQ':
            return(doc_field.value)
        if doc_field.VM == 1 and string:
            return(str(doc_field.value))
        else:
            return(doc_field.value)
    return(default)

def get_dataset(in_dir, case, suffix):
    filename = os.path.join(in_dir, case + suffix)
    img = None
    if os.path.exists(filename):
        img = sitk.ReadImage(filename)
        img = sitk.GetArrayFromImage(img)

    return img


def get_caselist(txt_file):
    datasets = []
    with open(txt_file, 'r') as f:
        for dset in f:
            datasets.append(dset.strip())
    return datasets

def read_dicom(in_dir, case, list_dcm):
    
# =============================================================================
#     dicom_names = os.listdir(os.path.join(in_dir, case))
#     dicom_paths = []
#     for dcm in dicom_names:
#         if dcm[:2] == 'CT':
#             dicom_paths.append(os.path.join(in_dir, case, dcm))
# =============================================================================

    dicom_paths = list_dcm
    img_positions = []
    for dcm in dicom_paths:
        ds = dcmread(dcm)
        img_positions.append(ds.ImagePositionPatient[2])

    indexes = np.argsort(np.asarray(img_positions))
    dicom_names = list(np.asarray(dicom_paths)[indexes])

    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(dicom_names)
    img = reader.Execute()

    return img


def read_dicom_dose(in_dir, case, list_dcm):
# =============================================================================
#     dicom_names = os.listdir(os.path.join(in_dir, case))
#     dicom_paths = []
#     for dcm in dicom_names:
#         if dcm[:2] == 'RD':
#             dose_file_name = os.path.join(in_dir, case, dcm)
# =============================================================================

    dose_file_name = list_dcm[0]
    #img_positions = []
    #for dcm in dicom_paths:
    #    ds = dcmread(dcm)
    # img = ds.pixel_array
    dose_img = dcmread(dose_file_name)
    # dose_img = sitk.ReadImage([dcm for dcm in dicom_paths])
    arr_dose = dose_img.pixel_array
    rt_dose = arr_dose*dose_img.DoseGridScaling
    rt_dose_itk = sitk.GetImageFromArray(rt_dose)
    rt_dose_itk.SetOrigin(dose_img.ImagePositionPatient)
    rt_dose_itk.SetSpacing([float(dose_img.PixelSpacing[0]), float(dose_img.PixelSpacing[1]), dose_img.GridFrameOffsetVector[1]-dose_img.GridFrameOffsetVector[0]])
    return rt_dose_itk


def get_rtstruct_dicom(in_dir, case, list_dcm):

# =============================================================================
#     rt_file = None
#     files = os.listdir(os.path.join(in_dir, case))
#     #pdb.set_trace()
#     for file in files:
#         if file[:2] == 'RS':
#             rt_file = file
# 
#     if rt_file is None:
#         return None
# =============================================================================

    rt_file = list_dcm
   ## rt_file = os.path.join(in_dir, case, rt_file)
    ds = dcmread(rt_file)

    return ds

def get_rtplan_dicom(in_dir, case, list_dcm):

# =============================================================================
#     rt_file = None
#     files = os.listdir(os.path.join(in_dir, case))
#     #pdb.set_trace()
#     for file in files:
#         if file[:2] == 'RS':
#             rt_file = file
# 
#     if rt_file is None:
#         return None
# =============================================================================

    rt_file = list_dcm
   ## rt_file = os.path.join(in_dir, case, rt_file)
    ds = dcmread(rt_file)
    return ds

def get_ref_ct(ct_in_dir, case, idx2):
    filename = os.path.join(ct_in_dir, case + '_'+str(idx2) +'_CT.nrrd')
    ref_ct = sitk.ReadImage(filename)

    return ref_ct


def get_oar_roi_indexes(ds, oar_dict):
    for idx, struc in enumerate(ds.StructureSetROISequence):
        elem = struc['ROIName']
        if elem.value in oar_dict:
            oar_dict[elem.value] = idx

    return oar_dict

def get_oar_roi_indexes_modified(ds, oar_dict, Plan_Dx):
    ##!! LY  important to find ROI names
    new_oar_dict = {}

## create a list of structure names in plan

    struc_name_plan = [];
    for idx, struc in enumerate(ds.StructureSetROISequence):
        elem = struc['ROIName']
        roi_name = elem.value.lower()
        struc_name_plan.append(roi_name)
           
    # Basic case ROI name difference is only upper or lower case
    for k, v in oar_dict.items():
        # We are searching for variations of key 'k' in the rt structure dicom data
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            if k == roi_name:
                new_oar_dict[k] = k
                oar_dict[k] = idx  # It means ROI index found for this anatomy
                break

    # Second basic case ROI name is for e.g. 'Cord1' or 'heart1' and similar variations with suffix '1'
    for k, v in oar_dict.items():
        # We are searching for variations of key 'k' in the rt structure dicom data
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            if roi_name == k+'1':
                new_oar_dict[k] = roi_name
                oar_dict[k] = idx  # It means ROI index found for this anatomy
                break


    # Third special case for Cord which can be sometimes named 'SpinalCord'
    v = oar_dict['lung_l']
    if v == -1:  # If still not found
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            if re.fullmatch(r'^l.*lung$|^lung.*l$!^lung.*left$', roi_name):
                new_oar_dict['lung_l'] = roi_name
                oar_dict['lung_l'] = idx  # It means ROI index found for this anatomy
                break               


    v = oar_dict['lung_r']
    if v == -1:  # If still not found
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            if re.fullmatch(r'^r.*lung$|^lung.*r$!^lung.*right$', roi_name):
                new_oar_dict['lung_r'] = roi_name
                oar_dict['lung_r'] = idx  # It means ROI index found for this anatomy
                break               

    v = oar_dict['spinalcord']
    if v == -1:  # If still not found
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            if re.fullmatch(r'cord', roi_name):
                new_oar_dict['spinalcord'] = roi_name
                oar_dict['spinalcord'] = idx  # It means ROI index found for this anatomy
            elif ('canal' in roi_name) & ('prv' not in roi_name):
                new_oar_dict['spinalcord'] = roi_name
                oar_dict['spinalcord'] = idx            
                break               
            else:
                pass
    v = oar_dict['ptv']
    if v == -1:  # If still not found
        for idx, struc in enumerate(ds.StructureSetROISequence):
            elem = struc['ROIName']
            roi_name = elem.value.lower()
            dx_str1 = "{:.0f}".format(Plan_Dx).rjust(2,'0')
            dx_str2 = "{:.0f}".format(Plan_Dx*100).rjust(4,'0')
            pat1 = re.compile(r'^ptv.?' + dx_str2 + '$|' + '^ptv.?' + dx_str2 + 'cGy$|' + '^ptv.?' + dx_str1 + '$|' \
                              + '^ptv.*eval$|^ptv.*final$|^ptv.*total$|^ptv.*all$|^ptv_\d+$|^ptv\d+$|^ptv.?\d+.?Gy$', re.IGNORECASE)
            if re.fullmatch(pat1, roi_name):
                new_oar_dict['ptv'] = roi_name
                oar_dict['ptv'] = idx  # It means ROI index found for this anatomy
                break               
            
    # Remaining cases: find any ROI names that starts with target anatomy name for e.g. 'PTV_New' 'PTV_Primary' etc.
    # Will analyze these cases manually
    for k, v in oar_dict.items():
        # We are searching for variations of key 'k' in the rt structure dicom data
        if v == -1:  # If still not found
            for idx, struc in enumerate(ds.StructureSetROISequence):
                elem = struc['ROIName']
                roi_name = elem.value.lower()
                if roi_name != k+'1':
                    if (roi_name in k) | (k in roi_name):
                        new_oar_dict[k] = str(idx)

    # Based on the analysis of the previous case output, if the length of the new dictionary is 6 then all anatomies
    # were found just some names were different. In this case just transfer the roi index to oar_dict
    # If the lenght of new dictionary is greater than 6 then multiple PTVs were found. Analysis of the dicom data in
    # slicer3d shows that in such cases we can pick ROI name 'E_AllPTVs' which combines all PTV masks into one.
    # If length of new dictionary is less than 6 then that is one special case '38068983' which has lungs named as
    # left_lung and right_lung. Find those id's and use them in oar_dict.
    # These above 3 cases based on length of dictionary should give us all RT Structs for all cases.

    # if len(new_oar_dict) == 6:
    #     for k, v in oar_dict.items():
    #         if v == -1:
    #             for new_k, new_v in new_oar_dict.items():
    #                 # if new_k.startswith(k):
    #                 #     oar_dict[k] = new_v
    #                 #     break
    #                 if k in new_k:
    #                     oar_dict[k] = new_v
    #                     break
    # elif len(new_oar_dict) > 6:
    #     for idx, struc in enumerate(ds.StructureSetROISequence):
    #         elem = struc['ROIName']
    #         roi_name = elem.value.lower()
    #         if roi_name == 'e_allptvs':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'ptv_total':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'ptv_60' or 'ptv60':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'ptvdvh':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'ptv_dibh6000':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'z_ptvoptr':
    #             oar_dict['ptv'] = idx
    #         elif roi_name == 'ptv final':
    #             oar_dict['ptv'] = idx

    # else:
    #     for idx, struc in enumerate(ds.StructureSetROISequence):
    #         elem = struc['ROIName']
    #         roi_name = elem.value.lower()
    #         if roi_name == 'left_lung':
    #             oar_dict['lung_l'] = idx
    #         if roi_name == 'right_lung':
    #             oar_dict['lung_r'] = idx
    #         if roi_name == '1ptv_rt_lung6000':
    #             oar_dict['ptv'] = idx

    return new_oar_dict, oar_dict


def get_contour_points(ds_elem):
    # Get the closed contour points triplets from dataset elem for the given plane
    points = np.asarray(ds_elem.ContourData)  # Array as [x0, y0, z0, x1, y1, z1 ....] in physical units
    points = np.reshape(points, (-1, 3))  # Array of point triplets, [ [x0, y0, z0], [x1, y1, z1], ....]

    return points


def transform_phy_pts_to_indexes(points, ref_ct):
    # Transform physical points to pixel coordinates (array indexes) using the reference ct's physical params
    pixel_coords = np.zeros_like(points)
    for idx in range(points.shape[0]):
        coord = (points[idx][0], points[idx][1], points[idx][2])
        coord = ref_ct.TransformPhysicalPointToIndex(coord)
        pixel_coords[idx] = coord

    return pixel_coords.astype(int)


def fill_planar_contour_as_mask(pixel_coords, mask3d):

    arr = np.zeros((mask3d.shape[1], mask3d.shape[2]), dtype='int32')  # 2D shape
    poly = pixel_coords[:,:2]
    #print(poly.dtype, arr.dtype)
    img = cv2.fillPoly(img=arr, pts=np.int32([poly]), color=1)  # Polygon fill the planar contour#Gourav changed it np.int32 as arr and points had different type
    mask = img.astype(np.uint8)
    zcord = pixel_coords[0][2]
    mask3d[zcord] = mask3d[zcord] + mask  # Single slice can have multiple contours which are specified separately (PTV mostly)

    return mask3d

def resample(img, ref_image):
    resampler = sitk.ResampleImageFilter()
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetReferenceImage(ref_image)
    img = resampler.Execute(img)
    
    return img

def write_image(img_arr, out_dir, case, suffix, ref_ct):
    img_itk = sitk.GetImageFromArray(img_arr)
    img_itk.SetOrigin(ref_ct.GetOrigin())
    img_itk.SetSpacing(ref_ct.GetSpacing())
    img_itk.SetDirection(ref_ct.GetDirection())
    filename = os.path.join(out_dir, case + suffix)
    sitk.WriteImage(img_itk, filename)

def get_crop_settings(oar, labels2):
    # Use to get crop settings
    # Don't use cord or eso as they spread through more slices
    # If total number of slices is less than 128 then don't crop at all
    # Use start and end index from presence of any anatomy or ptv
    # If that totals more than 128 slices then leave as is.
    # If that totals less than 128 slices then add slices before and after to make total slices to 128

    oar1 = oar.copy()
    oar1[np.where(oar == labels2['spinalcord'])] = 0
#    oar1[np.where(oar == labels2['spinalcord_prv'])] = 0
#    oar1[np.where(oar == labels2['esophagus_prv'])] = 0
    oar1[np.where(oar == labels2['esophagus'])] = 0
    oar1[np.where(oar == labels2['esophagus_ce'])] = 0

    # For 2D cropping just do center cropping 256x256
    center = [0, oar.shape[1] // 2, oar1.shape[2] // 2]
    start = [0, center[1] - 150, center[2] - 150]
    end = [0, center[1] + 150, center[2] + 150]

    depth = oar1.shape[0]
    if depth < 128:
        start[0] = 0
        end[0] = depth

        return start, end

    first_slice = -1
    last_slice = -1
    for i in range(depth):
        frame = oar1[i]
        if np.any(frame):
            first_slice = i
            break
    for i in range(depth-1, -1, -1):
        frame = oar1[i]
        if np.any(frame):
            last_slice = i
            break

    expanse = last_slice - first_slice + 1
    if expanse >= 128:
        start[0] = first_slice
        end[0] = last_slice

        return start, end

    #print('Get\'s here')
    slices_needed = 128 - expanse
    end_slices = slices_needed // 2
    beg_slices = slices_needed - end_slices

    room_available = depth - expanse
    end_room_available = depth - last_slice - 1
    beg_room_available = first_slice

    leftover_beg = beg_room_available - beg_slices
    if leftover_beg < 0:
        end_slices += np.abs(leftover_beg)
        first_slice = 0
    else:
        first_slice = first_slice - beg_slices

    leftover_end = end_room_available - end_slices
    if leftover_end < 0:
        first_slice -= np.abs(leftover_end)
        last_slice = depth - 1
    else:
        last_slice = last_slice + end_slices

    if first_slice < 0:
        first_slice = 0

    start[0] = first_slice
    end[0] = last_slice

    return start, end


def crop_resize_img(img, start, end, is_mask=False):
    # Crop to setting given by start/end coordinates list, assuming depth,height,width

    img_cropped = img[start[0]:end[0]+1, start[1]:end[1], start[2]:end[2]]
    img_cropped = np.moveaxis(img_cropped, 0, -1)  # Slices last

    order = 0
    if is_mask is False:
        order = 1
    img_resized = transform.resize(img_cropped, (128,128,128), order=order, preserve_range=True, anti_aliasing=False).astype(np.float32)
  ## not working for CV2 ignore  
  ##img_resized = cv2.resize(img_cropped, (128, 128, 128), interpolation=cv2.INTER_LINEAR).astype(
   #     np.float32)
    if is_mask is True:
        img_resized = img_resized.astype(np.uint8)

    img_resized = np.moveaxis(img_resized, -1, 0)  # Slices first again

    return img_resized


def get_torch_tensor(npy_tensor, device):
    out = torch.from_numpy(npy_tensor)
    out.to(device)

    return out

def get_dvh(dose, oar, ptv):
    # Compute and return the dvh for all 6 OAR structures
    #device = torch.device('cuda:0')
    device = torch.device('cpu')
    dose = get_torch_tensor(dose, device)
    oar = get_torch_tensor(oar, device).long()
    oar = torch.nn.functional.one_hot(oar, 8)[..., 1:]  # Remove BG
    #oar = oar.permute(3, 0, 1, 2).to(torch.float)
    oar = oar.permute(3, 0, 1, 2)
    ptv = get_torch_tensor(ptv, device).long().unsqueeze(dim=0)
    #ptv = ptv.to(torch.float)
    oar = torch.cat((oar, ptv), axis=0)

    vols = torch.sum(oar, axis=(1, 2, 3))
    n_bins = 376
    hist = torch.zeros((n_bins, 8)).to(device)
    bins = torch.linspace(0, 75, n_bins)
    bin_w = bins[1] - bins[0]
    
    ##return vols, bins
    for i in range(bins.shape[0]):
        diff = torch.sigmoid((dose - bins[i]) / bin_w)
        diff = torch.cat(8 * [diff.unsqueeze(axis=0)]) * oar
        num = torch.sum(diff, axis=(1, 2, 3))
        hist[i] = (num / vols)

## get DVH for Esophagus and A_LAD
    hist1 = torch.zeros((n_bins, 1)).to(device)
    hist2 = torch.zeros((n_bins, 1)).to(device)  
    oar_esoph =  oar[1,...]| oar[2,...]
    oar_heart =  oar[3,...]| oar[4,...]
    vol_esoph = torch.sum(oar_esoph)
    vol_heart = torch.sum(oar_heart)
    
    for i in range(bins.shape[0]):
        diff = torch.sigmoid((dose - bins[i]) / bin_w)
        diff1 = diff*oar_esoph
        num = torch.sum(diff1)
        hist1[i] = (num / vol_esoph)
        diff2 = diff*oar_heart
        num = torch.sum(diff2)
        hist2[i] = (num / vol_heart)
        
    #return hist, hist1
    hist.select_scatter(torch.squeeze(hist1),1,1)
    hist.select_scatter(torch.squeeze(hist2),1,3)
   # hist[:,1] = hist1
   # hist[:,3] = hist2

    hist_numpy = hist.cpu().numpy()
    bins_np = bins.cpu().numpy()

    return hist_numpy, bins_np

def get_dvh_clinical(dose, oar, ptv):
    # Compute and return the dvh for all 8 OAR structures
    # calculate real DVH using heavyside function instead of sigmoid
    
    device = torch.device('cpu')
    #device = torch.device('cuda:0')
    dose = get_torch_tensor(dose, device)
    oar = get_torch_tensor(oar, device).long()
    oar = torch.nn.functional.one_hot(oar, 8)[..., 1:]  # Remove BG
    oar = oar.permute(3, 0, 1, 2)
    ptv = get_torch_tensor(ptv, device).long().unsqueeze(dim=0)
    #ptv = ptv.to(torch.float)
    oar = torch.cat((oar, ptv), axis=0)

    vols = torch.sum(oar, axis=(1, 2, 3))
    n_bins = 376
    hist = torch.zeros((n_bins, 8)).to(device)
    bins = torch.linspace(0, 75, n_bins)
    bin_w = bins[1] - bins[0]

    for i in range(bins.shape[0]):
        diff = torch.heaviside((dose - bins[i]) / bin_w,torch.ones_like(dose))
        diff = torch.cat(8 * [diff.unsqueeze(axis=0)]) * oar
        num = torch.sum(diff, axis=(1, 2, 3))
        hist[i] = (num / vols)

## get DVH for Esophagus and A_LAD
    hist1 = torch.zeros((n_bins, 1)).to(device)
    hist2 = torch.zeros((n_bins, 1)).to(device)  
    oar_esoph =  oar[1,...]| oar[2,...]
    oar_heart =  oar[3,...]| oar[4,...]
    vol_esoph = torch.sum(oar_esoph)
    vol_heart = torch.sum(oar_heart)
    
    for i in range(bins.shape[0]):
        diff = torch.heaviside((dose - bins[i]) / bin_w, torch.ones_like(dose))
        diff1 = diff*oar_esoph
        num = torch.sum(diff1)
        hist1[i] = (num / vol_esoph)
        diff2 = diff*oar_heart
        num = torch.sum(diff2)
        hist2[i] = (num / vol_heart)
        
    hist.select_scatter(torch.squeeze(hist1),1,1)
    hist.select_scatter(torch.squeeze(hist2),1,3)
    #hist[:,1] = hist1
    #hist[:,3] = hist2

    hist_numpy = hist.cpu().numpy()
    bins_np = bins.cpu().numpy()

    return hist_numpy, bins_np
    
def process_case(in_dir, out_dir, case, idx2, labels2):
    ##!! LY importanct  Process case by case from each subfolder
    ct = get_dataset(out_dir, case+'_'+str(idx2), '_CT.nrrd')
    dose = get_dataset(out_dir, case+'_'+str(idx2), '_dose_resampled.nrrd')
    # dose = get_dataset(out_dir, case, '_dose_resampled.nrrd')  # Manual dose
    oar = get_dataset(out_dir, case+'_'+str(idx2), '_RTSTRUCTS.nrrd')
    ptv = get_dataset(out_dir, case+'_'+str(idx2), '_PTV.nrrd')
    # beamlet = get_dataset(out_dir, case, '_echo_dose_beamlet_resampled.nrrd')
    # beamlet = get_dataset(out_dir, case, '_manual_dose_beamlet_resampled.nrrd')
    # beamlet = get_dataset(out_dir, case, '_echo_dose_beamlet_sparse_resampled.nrrd')
    # beamlet = get_dataset(out_dir, case, '_manual_dose_beamlet_sparse_resampled.nrrd')

    oar_copy = oar.copy()
    oar_copy[np.where(ptv == 1)] = labels2['ptv']

    start, end = get_crop_settings(oar_copy, labels2)

    #print(start,end)
    ct = crop_resize_img(ct, start, end, is_mask=False)
    oar = crop_resize_img(oar, start, end, is_mask=True)
    ptv = crop_resize_img(ptv, start, end, is_mask=True)
    dose = crop_resize_img(dose, start, end, is_mask=False)
    # beamlet = crop_resize_img(beamlet, start, end, is_mask=False)
    # beamlet[np.where(ptv == 1)] = 60  # PTV volume set to prescribed dose

    # Scale PTV volume (region) in dose to have average prescibed 60 Gy

    #print(ct)   ## debug
    num_ptv = np.sum(ptv)
    dose_copy = dose.copy()
    dose_copy *= ptv
    #return ptv, oar
    #print(ptv)
    #print(dose_copy)
    sum = np.sum(dose_copy)
    scale_factor = (60 * num_ptv) / sum

    dose_copy *= scale_factor

    dose[np.where(ptv == 1)] = dose_copy[np.where(ptv == 1)]

    ct = np.clip(ct, a_min=-1000, a_max=3071)
    ct = (ct + 1000) / 4071
    ct = ct.astype(np.float32)

    dose = np.clip(dose, a_min=0, a_max=75)
    ##return oar, ptv

    hist_sig, bins = get_dvh(dose, oar, ptv)
    hist_clinical, bins = get_dvh_clinical(dose, oar, ptv)

    filename = os.path.join(out_dir, case+'_'+str(idx2) )
     
    return hist_clinical, bins
    

def batch_process_cases(in_dir, out_dir):
    # Initialize parser
# =============================================================================
#     parser = argparse.ArgumentParser()
#     
#     # Adding optional argument
#     parser.add_argument("--in_dir", required=True, help="Enter input dir having patient folders with their dicoms")
#     parser.add_argument("--out_dir", required=True, help="Enter out dir having patient folders with their dicoms")
#     args, _ = parser.parse_known_args()
#     in_dir = args.in_dir
#     out_dir = args.out_dir
# =============================================================================
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cases = os.listdir(in_dir)
    
    labels = {
        'spinalcord': 1,
        'esophagus': 2,
        'heart': 3,
        'lung_l': 4,
        'lung_r': 5,
        'ptv': 1
    }  # PTV will be stored separately as its extent is not mutually exclusive with other anatomies
    

            #oars = ['Cord', 'Esophagus', 'Heart', 'Lung_L', 'Lung_R', 'PTV']
          ##  oars = ['lungs','lung_l', 'lung_r','spinalcord', 'spinalcord_prv', 'esophagus', 'esophagus_prv', 'esophagus_ce','heart', 'a_lad',  'ptv']
            ## for 2024 AAPM ignore PRVs for now

    oars = ['lungs','lung_l', 'lung_r','spinalcord', 'esophagus', 'esophagus_ce','heart', 'a_lad',  'ptv']

    oars2 = ['lung_l', 'lung_r','spinalcord', 'esophagus', 'heart',  'ptv']    # for DL
    labels2 = dict.fromkeys(oars,0)
    label_seq = [0, 6, 7, 1, 2, 3, 4, 5, 1]
    for w2, key2 in enumerate(oars):
        labels2[key2] = label_seq[w2]        ## labels2 generate full organ mask to caluclate evalution DVHs not for learning  
        
    Labels3 = {key: labels2[key] for key in oars2}    ## lables3 for DL
    
## Data frame of case info
    column_names = ['CaseID', 'PlanInd', 'PlanID', 'Dx','AllDataSet?','AllContours?','ExtractedOK?'] + oars
    df_cases = pd.DataFrame(columns=column_names)
    filename_caseinfo = os.path.join(out_dir, 'Cases_processing_summary.csv')
    #fh = open(filename,"w",newline="")

    #for idx, case in enumerate(cases):
    for idx in range(0, len(cases)-1):
        case = cases[idx]
        if not case.endswith('.npz'):
            continue
        print('Processing case {}: {} of {} ...'.format(case, idx+1, len(cases)))
        case_path = os.path.join(in_dir, case)
        caseID = re.search(r'^VCU_Lung_\d+_\d', case).group()
        row_data1 = {'CaseID':caseID, 'PlanInd':0, 'PlanID': 0, 'Dx':0,'AllDataSet?': True}
        try:
            dict_np = np.load(case_path)
            w11 = dict_np['OAR']
            ptv = dict_np['PTV']
            w12 = w11.copy()
            if np.sum(np.where(w11==3))> 0.5:
                w12[np.where(w11==3)]=2
            if np.sum(np.where(w11==4))> 0.5:    
                w12[np.where(w11==4)]=3
            if np.sum(np.where(w11==5))> 0.5:
                w12[np.where(w11==5)]=3
            if np.sum(np.where(w11==6))> 0.5:
                w12[np.where(w11==6)]=4
            if np.sum(np.where(w11==7))> 0.5:
                w12[np.where(w11==7)]=5
            
            w13 = dict_np['HIST_SIG']
            w14 = np.delete(w13,[3,5],axis=1)
            w15 = dict_np['HIST_CLINICAL']
            w16 = np.delete(w15,[3,5],axis = 1)
            
            device = torch.device('cpu')
            #device = torch.device('cuda:0')
#=============================================================================
            w17 = get_torch_tensor(w12, device).long()
            w18 = torch.nn.functional.one_hot(w17, 6)
            w18 = w18.permute(3, 0, 1, 2)
            #!!ptv = get_torch_tensor(ptv, device).long().unsqueeze(dim=0)
            #!!oar_new = torch.cat((w18, ptv), axis=0)
#=============================================================================
            
            
            w19 = torch.sum(w18[1:,...], axis=(1, 2, 3))
            
            vols_red = dict.fromkeys(oars2[:-1],0)
            for i6, oar_name in enumerate(oars2):
                vols_red[oar_name] = w19[labels[oar_name]-1].item()
                
            row_data = {**row_data1, **vols_red}
            df_cases = df_cases.append(row_data,ignore_index=True)   
            
            if (vols_red['lung_l'] ==0) or (vols_red['lung_r'] == 0):
                filename_npz = os.path.join(out_dir, 'DNU_'+caseID+'_reduced.npz')
            else:
                filename_npz = os.path.join(out_dir, caseID+'_reduced.npz')

            np.savez(filename_npz, CT=dict_np['CT'], DOSE=dict_np['DOSE'], OAR=w12, PTV=dict_np['PTV'], HIST_SIG=w14, HIST_CLINICAL = w16, BINS=dict_np['BINS'])
     
            # Export the DataFrame as a CSV spreadsheet
            fh = open(filename_caseinfo,"w",newline="")
            df_cases.to_csv(fh)
            fh.flush()
            fh.close()
        except:
            print('Processing case {} failed...'.format(caseID))
    #fh.close()    
    # debug
    return df_cases, w16
## DEbug
#out_dir = "Z:\LulinY\Lung-dosimetrics\2024\python_code\DoseRTX\VCU_Lung_2024_processed_data"
in_dir = "C:\Lulin-home\KBP-lung\CE project\DoseRTX\VCU_Lung_2024_processed_data"
out_dir = "C:\\Lulin-home\\KBP-lung\\CE project\DoseRTX\\VCU_Lung_2024_dataset\\test"
w2, hist = batch_process_cases(in_dir,out_dir)
#fh.close()