#!/usr/bin/env python3
import pydicom
#import vtk
import numpy as np
#from scipy import ndimage
#from matplotlib.path import Path
#import vtk.util.numpy_support as vtk_np
#from bisect import bisect_left
import warnings
import json
import copy
import logging
import os

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)#DEBUG)
logger = logging.getLogger(__name__)

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

def dicom_validation(filepaths, recursive=True, strict=True, patient_id=None):
    """
    Organizes and validates DICOM-RT files into coherent groupings for radiotherapy data processing.
    
    This function processes a collection of potential DICOM files, identifying and grouping related 
    RTPLAN, RTDOSE, RTSTRUCT, RTRECORD, and RTIMAGE files based on cross-referential tags within the 
    DICOM metadata. Specifically, it ensures that each grouping contains a coherent set of files for 
    a single radiotherapy treatment plan for a single patient. The function performs preliminary validation 
    to ensure that files are genuine DICOM files, and then further checks for the presence of essential 
    components like RTPLAN, RTDOSE, and RTSTRUCT in each grouping. If a specified patient ID is provided, 
    the function will only process files that match the given ID, providing a way to filter the input 
    and ensure data consistency. Invalid, duplicate, or unmatched files, or those belonging to other/multiple 
    patients, are flagged with warnings, ensuring robustness in handling diverse inputs. Note that filepaths 
    are required/assumed to have been pre-assessed via security scanning software and have been deemed 'safe'.
    
    Args:
        filepaths (list of str): A list of potential DICOM file paths to be organized and validated.
        recursive (boolean, optional): An optional boolean to determine how to navigate directory input(s)
                                              If True (default), will include directory contents in processing.
                                              If False, will omit directory
        patient_id (str, optional): An optional patient ID to filter the processed files. 
                                              If provided, only files matching this ID will be considered.

    Returns:
        list of dict: A list of dictionaries, with each dictionary representing a validated grouping of 
                      related DICOM files. Each grouping contains references to the relevant DICOM files 
                      for RTPLAN, RTDOSE, RTSTRUCT, and optionally one or more RTRECORDs, one or more RTIMAGEs, and one or more IMAGEs.
        str: The patient ID associated with the processed DICOM files. If multiple patient IDs are 
             detected, the function will raise an exception if patient_id not prespecified.
        list of dict: A list of dictionaries corresponding to excluded files, with each dictionary representing a file-specific error. 
                      Each dictionary contains the file path and error type.

    Example usage:
        groupings, patient_id, errors = dicom_validation(filepaths, patient_id)
    """
    groupings = []
    excluded_files = []
    # dictionaries mapping IDs to files, and outgoing DICOM object reference(s)
    RTplans = {}
    RTdoses = {}
    RTstructs = {}
    RTrecords = {}
    RTimages = {}
    # we also want to identify image files associated with DICOM-RT plans - i.e. used for structure set / registration / etc. --> can link that then to the patient plan
    images = {}
    # TO DO: --> Get Tianjun GammaKnife RT Record data
    # TO DO:  get Reflexion data RT RECORD for Tianjun?? @ UT SW
    # TO DO: Lulin to get a Brainlab plan as well for comparison
    # TO DO: Rishabh to get a whole suite of test data from IHERO
    sop_classes = {
        '1.2.840.10008.5.1.4.1.1.481.1': 'RTIMAGE', # RT Image Storage
        '1.2.840.10008.5.1.4.1.1.481.2': 'RTDOSE', # RT Dose Storage
        '1.2.840.10008.5.1.4.1.1.481.3': 'RTSTRUCT', # RT Structure Set Storage
        '1.2.840.10008.5.1.4.1.1.481.4': 'RTRECORD', # RT Beams Treatment Record Storage
        '1.2.840.10008.5.1.4.1.1.481.5': 'RTPLAN', # RT Plan Storage
        '1.2.840.10008.5.1.4.1.1.481.6': 'RTRECORD', # RT Brachy Treatment Record Storage
        '1.2.840.10008.5.1.4.1.1.481.7': 'RTRECORD', # RT Treatment Summary Record Storage
        '1.2.840.10008.5.1.4.1.1.481.8': 'RTPLAN', # RT Ion Plan Storage
        '1.2.840.10008.5.1.4.1.1.481.9': 'RTRECORD', # RT Ion Beams Treatment Record Storage
        '1.2.840.10008.5.1.4.1.1.481.18': 'RTRECORD', # Tomotherapeutic Radiation Record Storage
        '1.2.840.10008.5.1.4.1.1.481.20': 'RTRECORD', # Robotic Radiation Record Storage (e.g. CyberKnife)
        '1.2.840.10008.5.1.4.1.1.481.23': 'RTIMAGE', # Enhanced RT Image Storage
        '1.2.840.10008.5.1.4.1.1.128': 'PET', # Positron Emission Tomography Image Storage
        '1.2.840.10008.5.1.4.1.1.128.1': 'PET', # Legacy Converted Enhanced PET Image Storage
        '1.2.840.10008.5.1.4.1.1.130': 'PET', # Enhanced PET Image Storage
        '1.2.840.10008.5.1.4.1.1.2': 'CT', # CT Image Storage
        '1.2.840.10008.5.1.4.1.1.2.1': 'CT', # Enhanced CT Image Storage
        '1.2.840.10008.5.1.4.1.1.2.2': 'CT', # Legacy Converted Enhanced CT Image Storage
        '1.2.840.10008.5.1.4.1.1.4': 'MRI', # MR Image Storage
        '1.2.840.10008.5.1.4.1.1.4.1': 'MRI', # Enhanced MR Image Storage
        '1.2.840.10008.5.1.4.1.1.4.2': 'MRI', # MR Spectroscopy Storage
        '1.2.840.10008.5.1.4.1.1.4.3': 'MRI', # Enhanced MR Color Image Storage
        '1.2.840.10008.5.1.4.1.1.4.4': 'MRI', # Legacy Converted Enhanced MR Image Storage
        '1.2.840.10008.5.1.4.1.1.3.1': 'US', # Ultrasound Multi-frame Image Storage
        '1.2.840.10008.5.1.4.1.1.20': 'NM', # Nuclear Medicine Image Storage
        '1.2.840.10008.5.1.4.1.1.66.1': 'REGISTER', # Spatial Registration Storage
        '1.2.840.10008.5.1.4.1.1.66.3': 'REGISTER', # Deformable Spatial Registration Storage
        '1.2.246.352.70.1.70': 'RTPLAN', # RT Plan Varian 1 Storage (Ethos radiotherapy system)
        '1.2.246.352.70.1.71': 'RTRECORD' # RT Treatment Record Varian 1 Storage (Ethos radiotherapy system)
    }
    ptID = None
    if isinstance(filepaths, str):
        filepaths = [filepaths]
    for filepath in filepaths:
        filepath = os.path.normpath(filepath)
        # Handle directory input
        if os.path.isdir(filepath):
            # Include file contents recursively, if desired (recursive == True)
            if recursive:
                for root, _, files in os.walk(filepath):
                    for file in files:
                        new_filepath = os.path.join(root, file)
                        if new_filepath not in filepaths:
                            filepaths.append(new_filepath)
            # Exclude the directory itself
            excluded_files.append({'filepath': filepath, 'type': 'IsDirectory'})
            continue
        # Exclude non-DICOM file(s)
        if not pydicom.misc.is_dicom(filepath):
            excluded_files.append({'filepath': filepath, 'type': 'NonDICOM'})
            continue
        # Read DICOM metadata without loading pixel data
        try:
            dcm = pydicom.dcmread(filepath, stop_before_pixels=True)
        except:
            # Read headerless DICOM metadata without loading pixel data
            try:
                dcm = pydicom.dcmread(filepath, force=True, stop_before_pixels=True)
            # Exclude file(s) if cannot read DICOM metadata
            except:
                excluded_files.append({'filepath': filepath, 'type': 'InvalidDICOM'})
                continue
            else:
                # Exclude file(s) if SOPClassUID not present in DICOM metadata
                #(0008,0016) = SOP Class UID [Required (1), unique identifier]
                if 'SOPClassUID' not in dcm:
                    excluded_files.append({'filepath': filepath, 'type': 'InvalidDICOM'})
                    continue
        # Only process RTDOSE, RTPLAN, RTSTRUCT, RTRECORD, RTIMAGE file(s) + associated image(s) [CT, MR, PET, NM, US]
        #(0008,0016) = SOP Class UID [Required (1), unique identifier]
        sop_class = _get_tag(dcm,0x0008,0x0016)
        if sop_class not in sop_classes:
            excluded_files.append({'filepath': filepath, 'type': 'NonRTModule'})
            continue
        dcm_type = sop_classes[sop_class]
        if dcm_type == 'RTPLAN':
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            plan_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate RTPLAN(s)
            if plan_ID in RTplans:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateRTPLAN'})
                continue
            # Exclude rejected RTPLAN(s)
            #(300E,0002) = Approval Status [Required (1), code string]
            if strict and _get_tag(dcm,0x300E,0x0002) == 'REJECTED':
                excluded_files.append({'filepath': filepath, 'type': 'RejectedRTPLAN'})
                continue
            # Extract and reconcile patient identity (RTPLAN only - other files linked to RTPLAN)
            pat_ID = match_pt_identity_from_dicom(dcm)
            # Exclude file if patient ID anonymized
            if strict and pat_ID == 'ANONYMOUS':
                excluded_files.append({'filepath': filepath, 'type': 'AnonymizedRTPLAN'})
                continue
            if ptID is None:
                ptID = pat_ID
            # Exclude file(s) from other patient(s)
            if patient_id and patient_id != pat_ID:
                excluded_files.append({'filepath': filepath, 'type': 'WrongPatientRTPLAN'})
                continue
            # If RTPLANS for multiple patients without specifying desired patient, raise an exception
            elif ptID != pat_ID:
                pass
                #raise Exception(f'Multiple patients detected in the provided RTPLAN files. Mismatched IDs: {ptID} and {pat_ID}')
            #(300A,000C) = Plan Geometry [Required (1), code string, values:('PATIENT','TREATMENT_DEVICE')]
            RTSTRUCT_expected = (_get_tag(dcm,0x300A,0x000C,'PATIENT') == 'PATIENT')
            structset_IDs = []
            if RTSTRUCT_expected:
                structset_IDs.extend([
                    #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                    _get_tag(structset_ref,0x0008,0x1155)
                    #(300C,0060) = Referenced Structure Set Sequence [Conditional (1C), sequence]
                    for structset_ref in _get_tag(dcm,0x300C,0x0060,[])
                ])
            dose_IDs = []
            dose_IDs.extend([
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                _get_tag(dose_ref,0x0008,0x1155)
                #(300C,0080) = Referenced Dose Sequence [Optional (3), sequence]
                for dose_ref in _get_tag(dcm,0x300C,0x0080,[])
            ])
            #(0020,0052) = Frame of Reference UID [Required (1), unique identifier]
            frame_of_ref =  _get_tag(dcm,0x0020,0x0052)
            RTplans[plan_ID] = {'filepath': filepath, 'RTSTRUCT': structset_IDs, 'RTDOSE': dose_IDs, 'fref': frame_of_ref}
        elif dcm_type == 'RTDOSE':
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            dose_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate RTDOSE(s)
            if dose_ID in RTdoses:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateRTDOSE'})
                continue
            # Exclude RTDOSE(s) lacking dose grid
            #(7FE0,0010) = Pixel Data [Conditional (1C), OB or OW bit string]
            with open(filepath, 'rb') as f:
                if (0x7FE0, 0x0010) not in pydicom.filereader.read_partial(f, specific_tags=[(0x7FE0, 0x0010)], defer_size=1):
                    excluded_files.append({'filepath': filepath, 'type': 'RTDOSEMissingDoseGrid'})
                    continue
            # Exclude RTDOSE(s) lacking physical or effective dose information
            #(3004,0004) = Dose Type [Required (1), code string]
            if _get_tag(dcm,0x3004,0x0004) == 'ERROR':
                excluded_files.append({'filepath': filepath, 'type': 'RTDOSEMissingPhysOrEffDose'})
                continue
            # Exclude RTDOSE(s) encoded as relative doses
            #(3004,0002) = Dose Units [Required (1), code string]
            if strict and _get_tag(dcm,0x3004,0x0002) == 'RELATIVE':
                excluded_files.append({'filepath': filepath, 'type': 'RelativeNotAbsoluteRTDOSE'})
                continue
            #(3004,000A) = Dose Summation Type [Required (1), code string]
            dose_sum_type = _get_tag(dcm,0x3004,0x000A, "")
            plan_IDs = []
            record_IDs = []
            if dose_sum_type in ['PLAN','MULTI_PLAN','FRACTION','BEAM','BRACHY','FRACTION_SESSION','BEAM_SESSION','BRACHY_SESSION','CONTROL_POINT']:
                #TO DO!!!  FIGURE OUT re: BrainLab plan if de-linked from dose because of anonymization or is that different format encoding?!
                # TO DO!!!  convert to checking references via the (0008,1155) tag?!?!
                #temp = []
                #for series in _get_tag(dcm,0x0008,0x1115,[]):
                #    for instance in _get_tag(series,0x0008,0x114A,[]):
                #        if sop_classes[_get_tag(instance,0x0008,0x1150,'')] == 'RTPLAN':
                #            temp.extend([_get_tag(instance,0x0008,0x1155)]) 
                plan_IDs.extend([
                    #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                    _get_tag(refplan,0x0008,0x1155)
                    #(300C,0002) = Referenced RT Plan Sequence [Conditional (1C), sequence]
                    for refplan in _get_tag(dcm,0x300C,0x0002,[])
                ])
            elif dose_sum_type == 'RECORD':
                record_IDs.extend([
                    #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                    _get_tag(record, 0x0008, 0x1155)
                    #(3008,0030) = Referenced Treatment Record Sequence [Conditional (1C), sequence]
                    for record in _get_tag(dcm,0x3008,0x0030,[])
                ])
            else:
                excluded_files.append({'filepath': filepath, 'type': 'InvalidRTDOSE'})
                continue
            RTdoses[dose_ID] = {'filepath': filepath, 'RTPLAN': plan_IDs, 'RTRECORD': record_IDs}
        elif dcm_type == 'RTSTRUCT':
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            structset_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate RTSTRUCT(s)
            if structset_ID in RTstructs:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateRTSTRUCT'})
                continue
            #(0020,0052) = Frame of Reference UID [Required (1), unique identifier]
            frame_of_ref = _get_tag(dcm,0x0020,0x0052)
            RTstructs[structset_ID] = {'filepath': filepath, 'fref': frame_of_ref, 'IMAGES': []}
            #(3006,0010) = Referenced Frame of Reference Sequence [IHE Required (1), sequence]
            for rfref in _get_tag(dcm,0x3006,0x0010,[]):
                #(3006,0012) = RT Referenced Study Sequence [IHE Required (1), sequence]
                for study in _get_tag(rfref,0x3006,0x0012,[]):
                    #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]                
                    study_ID = _get_tag(study,0x0008,0x1155)
                    #(3006,0014) = RT Referenced Series Sequence [Required (1), sequence]
                    for series in _get_tag(study,0x3006,0x0014,[]):
                        #(0020,000E) = Series Instance UID [Required (1), unique identifier]                
                        series_ID = _get_tag(series,0x0020,0x000E)
                        #(3006,0016) = Contour Image Sequence [Required (1), sequence]
                        #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]                
                        image_IDs = [_get_tag(i,0x0008,0x1155) for i in _get_tag(series,0x3006,0x0016,[])]
                        RTstructs[structset_ID]['IMAGES'].append({'studyID': study_ID, 'seriesID': series_ID, 'images': image_IDs})
        elif dcm_type == 'RTRECORD':
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            record_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate RTRECORD(s)
            if record_ID in RTrecords:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateRTRECORD'})
                continue
            plan_IDs = []
            plan_IDs.extend([
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                _get_tag(refplan,0x0008,0x1155)
                #(300C,0002) = Referenced RT Plan Sequence [Required (2), sequence]
                for refplan in _get_tag(dcm,0x300C,0x0002,[])
            ])
            plan_IDs.extend([
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                _get_tag(instance,0x0008,0x1155)
                #(0008,1115) = Referenced Series Sequence [Conditional (1C), sequence]
                for series in _get_tag(dcm,0x0008,0x1115,[])
                #(0008,1115) = Referenced Instance Sequence [Required (1), sequence]
                for instance in _get_tag(series,0x0008,0x114A,[])
                #(0008,1150) = Referenced SOP Class UID [Required (1), unique identifier]
                if (
                    _get_tag(instance,0x0008,0x1150,'') in sop_classes and 
                    sop_classes[_get_tag(instance,0x0008,0x1150,'')] == 'RTPLAN'
                )
            ])
            record_IDs = []
            record_IDs.extend([
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                _get_tag(record,0x0008,0x1155)
                #(3008,0030) = Referenced Treatment Record Sequence [Optional (3), sequence]
                for record in _get_tag(dcm,0x3008,0x0030,[])
            ])
            # PERFORM SOME VALIDATION ON DOSE FILE
            RTrecords[record_ID] = {'filepath': filepath, 'RTPLAN': plan_IDs, 'RTRECORD': record_IDs}
        elif dcm_type == 'RTIMAGE':
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            image_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate RTIMAGE(s)
            if image_ID in RTimages:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateRTIMAGE'})
                continue
            refplan_ID = None
            #(300C,0002) = Referenced RT Plan Sequence [Optional (3), sequence]
            for refplan in _get_tag(dcm,0x300C,0x0002, []):
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                refplan_ID = _get_tag(refplan,0x0008,0x1155,None)
            RTimages[image_ID] = {'filepath': filepath, 'RTPLAN': refplan_ID}
        elif dcm_type in ['CT','MR','PET','NM','US']:
            # TO DO: check images to ensure that all slices are present and accounted for? - no easy way to do this using DICOM standard . . . could pull the instance number (0020,0013) and ensure that the range is complete & monotonically increasing with each new slice (if so, it's complete image) - but if not, doesn't mean it's NOT a complete image!
            #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
            image_ID = _get_tag(dcm,0x0008,0x0018,'')
            # Exclude duplicate CT/MR/PET/NM/US image(s)
            if image_ID in images:
                excluded_files.append({'filepath': filepath, 'type': 'DuplicateIMAGE'})
                continue
            #(0020,000D) = Study Instance UID [Required (1), unique identifier]
            study_ID = _get_tag(dcm,0x0020,0x000D)
            #(0020,000E) = Series Instance UID [Required (1), unique identifier]
            series_ID = _get_tag(dcm,0x0020,0x000E)
            #(0020,0052) = Frame of Reference UID [Required (1), unique identifier]
            frame_of_ref = _get_tag(dcm,0x0020,0x0052)
            images[image_ID] = {'filepath': filepath, 'studyID': study_ID, 'series_ID': series_ID, 'fref': frame_of_ref}
        else: #dcm_type == 'REGISTER'
            # TO DO -- perform linkages to ensure REGISTER is stored and linked, for instance, to DOSE GRID and appropriate other image components
            excluded_files.append({'filepath': filepath, 'type': 'RTRegistrationFile'})
            continue
    included_files = []
    for plan_ID, plan_data in RTplans.items():
        group = {'RTPLAN': plan_data['filepath'], 'RTDOSE': [], 'RTRECORD': [], 'RTIMAGE': [], 'IMAGE': []}
        # Exclude RTPLAN(s) without corresponding RTDOSE file(s)
        if any(d not in RTdoses for d in plan_data['RTDOSE']):
            excluded_files.append({'filepath': plan_data['filepath'], 'type': 'RTPLANMissingRTDOSE'})
            continue
        if not any(pid == plan_ID for pid in [item for d,dd in RTdoses.items() for item in dd['RTPLAN']]):
            excluded_files.append({'filepath': plan_data['filepath'], 'type': 'RTPLANMissingRTDOSE'})
            continue
        # Exclude RTPLAN(s) without corresponding RTSTRUCT file(s) if expected
        if any(s not in RTstructs for s in plan_data['RTSTRUCT']):
            excluded_files.append({'filepath': plan_data['filepath'], 'type': 'RTPLANMissingRTSTRUCT'})
            continue
        elif plan_data['RTSTRUCT']:
            # Assume just ONE RT structure set, but not verifying here
            struct_ID = plan_data['RTSTRUCT'][0]
            group['RTSTRUCT'] = RTstructs[struct_ID]['filepath']
            included_files.append(struct_ID)
            # Associate CT/MR/PET image(s) with RTSTRUCT via Frame of Reference UID
            if not RTstructs[struct_ID]['IMAGES']:
                for image_ID, image_data in images.items():
                    if image_data['fref'] == RTstructs[struct_ID]['fref']:
                        group['IMAGE'].append(image_data['filepath'])
                        included_files.append(image_ID)
            # Associate CT/MR/PET/NM/US image(s) with RTSTRUCT via linked image UIDs
            else:
                for i in [i['images'] for i in RTstructs[struct_ID]['IMAGES']]:  
                    # Exclude CT/MR/PET study if there are missing image(s)
                    if any(img not in images for img in i):
                        for img in [exclude for exclude in i if exclude in images]:
                            excluded_files.append({'filepath': images[img]['filepath'], 'type': 'IncompleteImage'})
                            del images[img]
                    else:
                        for img in i:
                            group['IMAGE'].append(images[img]['filepath'])
                            included_files.append(img)
        else:
            group['RTSTRUCT'] = None
            # Associate CT/MR/PET image(s) with RTPLAN
            for image_ID, image_data in images.items():
                if image_data['fref'] == plan_data['fref']:
                    group['IMAGE'].append(image_data['filepath'])
                    included_files.append(image_ID)
        # Associate RTRECORD(s) with RTPLAN
        for record_ID, record_data in RTrecords.items():
            if plan_ID in record_data['RTPLAN']:
                group['RTRECORD'].append(record_data['filepath'])
                included_files.append(record_ID)
        # Associate RTDOSE(s) with RTPLAN(s)
        for dose_ID, filepath in [(d,dd['filepath']) for d,dd in RTdoses.items() for pid in dd['RTPLAN'] if pid == plan_ID]:
            group['RTDOSE'].append(filepath)
            included_files.append(dose_ID)
        # Associate RTIMAGE(s) with RTPLAN
        for image_ID, image_data in RTimages.items():
            if image_data['RTPLAN'] == plan_ID:
                group['RTIMAGE'].append(image_data['filepath'])
                included_files.append(image_ID)
        groupings.append(group)
    # Exclude RTSTRUCT(s) without corresponding RTPLAN file(s)
    for filepath in [sdat['filepath'] for s,sdat in RTstructs.items() if s not in included_files]:
        excluded_files.append({'filepath': filepath, 'type': 'RTSTRUCTMissingRTPLAN'})
    # Exclude RTDOSE(s) without corresponding RTPLAN file(s)
    for filepath in [ddat['filepath'] for d,ddat in RTdoses.items() if d not in included_files]:
        excluded_files.append({'filepath': filepath, 'type': 'RTDOSEMissingRTPLAN'})
    # Exclude RTRECORD(s) without corresponding RTPLAN file(s)
    for filepath in [rdat['filepath'] for r,rdat in RTrecords.items() if r not in included_files]:
        excluded_files.append({'filepath': filepath, 'type': 'RTRECORDMissingRTPLAN'})
    # Exclude RTIMAGE(s) without corresponding RTPLAN file(s)
    for filepath in [idat['filepath'] for i,idat in RTimages.items() if i not in included_files]:
        excluded_files.append({'filepath': filepath, 'type': 'RTIMAGEMissingRTPLAN'})
    # Exclude CT/MR/PET image(s) without corresponding RTPLAN/RTSTRUCT file(s)
    for filepath in [idat['filepath'] for i,idat in images.items() if i not in included_files]:
        excluded_files.append({'filepath': filepath, 'type': 'ImageMissingRTObject'})
    return groupings, ptID, excluded_files

def match_pt_identity_from_dicom(dcm, source_info={}, patient_id_api=None, fuzzy_threshold=0.9):
    """
    Reconciles patient identity information extracted from a DICOM file's patient module 
    by leveraging API accessing patient identity matching service from VA's MPI.

    This function uses the provided patient identity attributes (such as name, birthdate, 
    sex, address information, etc.) to query an API that provides a list of potential 
    patient matches. Using fuzzy matching algorithms, the function identifies the most 
    likely match based on the provided DICOM data. If the match score exceeds a specified 
    threshold, the function returns the discrete individual patient ID from the API. If 
    no match exceeds the threshold, the function raises an exception.

    Args:
        dcm (object): A pydicom DICOM object.
        source_info (dict): A dictionary containing relevant source information pertaining to the DICOM
                             May include keys such as 'referral_id', 'sta3n', 'state', etc.
        patient_id_api (object): An API object or function that can be used to query a list 
                                 of potential patient matches based on the provided DICOM data.
        fuzzy_threshold (float, optional): The threshold for fuzzy matching algorithms to 
                                           consider a potential patient match as the correct match. 
                                           Defaults to 0.9.

    Returns:
        str: The discrete individual patient ID from the API that corresponds to the 
             patient information provided.

    Raises:
        Exception: If no match exceeds the specified threshold, or if there are other 
                   issues with the provided patient attributes or API query.
    """
    #(0012,0062) = Patient Identity Removed [Optional (3), code string, values:('YES','NO')]
    #(0012,0063) = De-identification Method [Conditional (1C), long string]
    #(0012,0064) = De-identification Method Code Sequence [Conditional (1C), sequence]
    if (_get_tag(dcm,0x0012,0x0062,'NO') == 'YES' or
            _get_tag(dcm,0x0012,0x0063,None) or
            len(_get_tag(dcm,0x0012,0x0064,[])) > 0):
        return('ANONYMOUS')
    #(0010,0010) = Patient's Name [Required (2), person name]
    pt_name = _get_tag(dcm,0x0010,0x0010,'',string=False)
    if pt_name:
        given_name = pt_name.given_name
        middle_name = pt_name.middle_name
        family_name = pt_name.family_name
        prefix = pt_name.name_prefix
        suffix = pt_name.name_suffix
    else:
        given_name = None
        middle_name = None
        family_name = None
        prefix = None
        suffix = None
    #(0010,0020) = Patient ID [Required (2), long string]
    ID = str(_get_tag(dcm,0x0010,0x0020,''))
    #(0010,0030) = Patient's Birth Date [Required (2), date (YYYYMMDD is DICOM standard)]
    DOB = _get_tag(dcm,0x0010,0x0030)
    if DOB:
        DOB = DOB[:4]+'-'+DOB[4:6]+'-'+DOB[6:8]
    #(0010,0040) = Patient's Sex [Required (2), code string]
    sex = _get_tag(dcm,0x0010,0x0040)
    #(0010,0101) = Other Patient Names [Optional (3), person name]
    name_other = _get_tag(dcm,0x0010,0x0101,'',string=False)
    #(0010,1002) = Other Patient IDs Sequence [Optional (3), sequence]
    #(0010,0020) = Patient ID [Required (1), long string]
    ID_other = [_get_tag(i,0x0010,0x0020) for i in _get_tag(dcm,0x0010,0x1002,[])]
    #(0010,2160) = Ethnic Group [Optional (3), short string]
    ethnicity = _get_tag(dcm,0x0010,0x2160,'')
    #(0008,0020) = Study Date [Required (2), date (YYYYMMDD is DICOM standard)]
    study_date = _get_tag(dcm,0x0008,0x0020) 
    # ^^ potentially use study date to identify pt because they should have some sort of RT activity in system around same time as study date
    # Unclear value of responsible person data for identity matching at this time
    #(0010,2297) = Responsible Person [Optional (2C), person name]
    #(0010,2298) = Responsible Person Role [Required (1C), code string]
    # Reconcile patient identity and map to ___ (?PatientICN, EDIPI, SSN, other?)
    return "<placeholder-ID>"

def dicom_extract(id=None, RTplan=None, RTstruct=None, RTdose=[], RTrecord=[]):
    """
    Extracts structured radiotherapy data from related DICOM-RT objects.

    This function processes the provided DICOM-RT files corresponding to a single plan 
    for a single patient. It retrieves relevant structured data, including dose-volume histograms, 
    treatment fractions, beam data, and more. This structured data is formatted for storage in a 
    MongoDB database.

    The function is built to process DICOM-RT files compliant with the DICOM standard, 
    specifically the RTPLAN, RTDOSE, RTSTRUCT, and RTRECORD SOP Classes. For detailed information 
    on these classes and their attributes, refer to: https://dicom.innolitics.com/ciods

    Args:
        id (str, optional): A user-specified patient identifier. If provided, it overrides the patient 
                            identifiers from the DICOM data.
        RTplan (str): File path to a single valid RTPLAN .dcm file. This file should be associated with 
                      the provided RTstruct and RTdose files. [REQUIRED]
        RTstruct (str): File path to a single valid RTSTRUCT .dcm file. This file should be associated 
                        with the provided RTplan and RTdose files. [REQUIRED (conditionally)]
        RTdose (list): A list of file paths to valid RTDOSE .dcm files associated with 
                      the provided RTplan and RTstruct files. [REQUIRED]
        RTrecord (list, optional): A list of file paths to valid RTRECORD .dcm files associated with the 
                                   provided RTplan file.

    Returns:
        str: A JSON-formatted string containing the extracted structured radiotherapy data. The structure 
             includes, but is not limited to:
             - Patient data (ID, name, birthdate)
             - Plan data (name, date, beams, modality)
             - Fraction data (fx groupings)
             - Dose data (max/min values, units)
             - Structure set (organs, ROIs, volumes)
             - DVH data (histograms, metrics)
             - Treatment records (sessions, control points, parameters)

    Note:
        This function presumes a set of already-validated and co-associated DICOM files, e.g. via dicom_validation()
    """
    if isinstance(RTdose, str):
        RTdose = [RTdose]
    if isinstance(RTrecord, str):
        RTrecord = [RTrecord]
    logger.info(f'Extracting DICOM-RT from RTPLAN: {RTplan}, RTSTRUCT: {RTstruct}, {len(RTdose)} RTDOSE file(s), and {len(RTrecord)} RTRECORD file(s)...')
    dcm = {key: {} for key in ['pt', 'plan', 'beams', 'fxs', 'structures', 
                                'ref_plans', 'treatments', 'brachy']}
    # Read in RTPLAN data
    plan_data = pydicom.dcmread(RTplan)
    # Read in RTSTUCT data (if expected)
    #(300A,000C) = Plan Geometry [Required (1), code string]
    if _get_tag(plan_data,0x300A,0x000C,'') == 'PATIENT':
        struct_data = pydicom.dcmread(RTstruct)
    else:
        logger.warning(f'RT Plan ({RTplan}) is not associated with a structure set')
        struct_data = None
    #------- RT Plan : Approval Module -------------------------------------------------------------
    #(300E,0002) = Approval Status [Required (1), code string]
    dcm['plan']['approval'] = _get_tag(plan_data,0x300E,0x0002)
    if dcm['plan']['approval'] != 'APPROVED':
        logger.warning(f'RT Plan ({RTplan}) is not an approved plan')
    #(300E,0004) = Review Date [Conditional (2C), date (YYYYMMDD is DICOM standard)]
    date = _get_tag(plan_data,0x300E,0x0004)
    if date:
        dcm['plan']['review_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
    else:
        dcm['plan']['review_date'] = None
    #(300E,0005) = Review Time [Conditional (2C), time (HHMMSS.FFFFFF&ZZXX is DICOM standard)]
    #(300E,0008) = Reviewer Name [Conditional (2C), person name]
    dcm['plan']['reviewer_name'] = _get_tag(plan_data,0x300E,0x0008)
    if (dcm['plan']['reviewer_name'] is not None):
        dcm['plan']['reviewer_name'] = str(dcm['plan']['reviewer_name'])
    #------- RT Plan : RT General Plan Module ------------------------------------------------------
    #(300A,0002) = RT Plan Label [Required (1), short string]
    dcm['plan']['plan_label'] = _get_tag(plan_data,0x300A,0x0002)
    #(300A,0003) = RT Plan Name [Optional (3), long string] 
    dcm['plan']['plan_name'] = _get_tag(plan_data,0x300A,0x0003) 
    #(300A,0004) = RT Plan Description [Optional (3), short text]
    dcm['plan']['plan_description'] = _get_tag(plan_data,0x300A,0x0004) 
    #(300A,0006) = RT Plan Date [Required (2), date (YYYYMMDD is DICOM standard)]
    date = _get_tag(plan_data,0x300A,0x0006)
    if date:
        dcm['plan']['plan_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
    else:
        dcm['plan']['plan_date'] = None
    #(300A,0007) = RT Plan Time [Required (2), time (HHMMSS.FFFFFF&ZZXX is DICOM standard)]
    time = (_get_tag(plan_data,0x300A,0x0007,'')[:6]+'000000')[:6]
    dcm['plan']['plan_time'] = time[:2]+':'+time[2:4]+':'+time[4:6]
    #(300A,000A) = Plan Intent [Optional (3), code string]
    dcm['plan']['plan_intent'] = _get_tag(plan_data,0x300A,0x000A)
    #(300A,000B) = Treatment Sites [Optional (3), long string]
    dcm['plan']['treatment_sites'] = _get_tag(plan_data,0x300A,0x000B)
    #------- RT Plan : General Equipment Module ----------------------------------------------------
    #(0008,0080) = Institution Name [Optional (3), long string]
    dcm['plan']['facility'] = _get_tag(plan_data,0x0008,0x0080)
    #------- RT Plan : SOP Common Module -----------------------------------------------------------
    #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
    plan_ID = _get_tag(plan_data,0x0008,0x0018)
    dcm['plan']['planID'] = plan_ID
    #------- RT Plan : General Study Module --------------------------------------------------------
    #(0008,0020) = Study Date [Required (2), date (YYYYMMDD is DICOM standard)]
    date = _get_tag(plan_data,0x0008,0x0020)
    if date:
        dcm['pt']['study_date'] = date[:4]+"-"+date[4:6]+"-"+date[6:8]
    else:
        dcm['pt']['study_date'] = None
    #(0008,1030) = Study Description [Optional (3), long string]
    dcm['pt']['study_description'] = _get_tag(plan_data,0x0008,0x1030)
    #(0008,1048) = Physician(s) of Record [Optional (3), person name]
    dcm['pt']['physicians'] = _get_tag(plan_data,0x0008,0x1048)
    if (dcm['pt']['physicians'] is not None):
        dcm['pt']['physicians'] = str(dcm['pt']['physicians'])
    #(0008,0090) = Referring Physician's Name [Required (2), person name]
    dcm['pt']['physician_referring'] = _get_tag(plan_data,0x0008,0x0090)
    if (dcm['pt']['physician_referring'] is not None):
        dcm['pt']['physician_referring'] = str(dcm['pt']['physician_referring'])
    #(0008,009C) = Consulting Physician's Name [Optional (3), person name]
    dcm['pt']['physician_consulting'] = _get_tag(plan_data,0x0008,0x009C)
    if (dcm['pt']['physician_consulting'] is not None):
        dcm['pt']['physician_consulting'] = str(dcm['pt']['physician_consulting'])
    #(0032,1033) = Requesting Service [Optional (3), long string]
    dcm['pt']['requesting_service'] = _get_tag(plan_data,0x0032,0x1033)
    #------- RT (Ion) Plan : RT (Ion) Beams Module -------------------------------------------------------------
    dcm['plan']['modality'] = None
    open_portfilms = 0
    trmt_portfilms = 0
    #(300A,00B0) = Beam Sequence [Required (1), sequence]
    #(300A,03A2) = Ion Beam Sequence [Required (1), sequence]
    for beam in list(_get_tag(plan_data,0x300A,0x00B0,[]))+list(_get_tag(plan_data,0x300A,0x03A2,[])):
        bmdata = {'planID': plan_ID, 'fractions': {}}
        #(300A,00C0) = Beam Number [Required (1), integer string]
        beam_ID = _get_tag(beam,0x300A,0x00C0)
        bmdata['beamID'] = beam_ID
        #(300C,006A) = Referenced Patient Setup Number [Optional (3), integer string]
        bmdata['ref_setupID'] = _get_tag(beam,0x300C,0x006A)
        #(300A,00C2) = Beam Name [Optional (3), long string]
        bmdata['beam_name'] = _get_tag(beam,0x300A,0x00C2)
        #(300A,00C3) = Beam Description [Optional (3), short text]
        bmdata['beam_description'] = _get_tag(beam,0x300A,0x00C3)
        #(300A,010E) = Final Cumulative Meterset Weight [Conditional (1C), decimal string]
        bmdata['cumulative_meterset'] = float(_get_tag(beam,0x300A,0x010E))
        if bmdata['cumulative_meterset'] == 100 or bmdata['cumulative_meterset'] == 1.0:
            bmdata['cumulative_meterset'] = None
        #(300A,00CE) = Treatment Delivery Type [Optional (3), code string]
        bmdata['delivery_type'] = _get_tag(beam,0x300A,0x00CE,'NONE')
        if bmdata['delivery_type'] not in ['TREATMENT','CONTINUATION']:
            if bmdata['delivery_type'] == 'OPEN_PORTFILM':
                open_portfilms += 1
            elif bmdata['delivery_type'] == 'TRMT_PORTFILM':
                trmt_portfilms += 1
            dcm['beams'][beam_ID] = bmdata
            continue
        #(300A,00C4) = Beam Type [Required (1), code string]
        bmdata['motion'] = _get_tag(beam,0x300A,0x00C4)
        #(300A,00C6) = Radiation Type [Required (2), code string]
        bmdata['modality'] = _get_tag(beam,0x300A,0x00C6)
        if not dcm['plan']['modality']:
            dcm['plan']['modality'] = bmdata['modality']
        elif dcm['plan']['modality'] != bmdata['modality']:
            dcm['plan']['modality'] = 'MIXED'
        #(300A,00C7) = High-Dose Technique Type [Conditional (1C), code string]
        bmdata['high_dose'] = _get_tag(beam,0x300A,0x00C7)
        #(300A,00B3) = Primary Dosimeter Unit [Optional (3), code string]
        bmdata['dosimeter_units'] = _get_tag(beam,0x300A,0x00B3)
        #(300A,00D0) = Number of Wedges [Required (1), integer string]
        bmdata['wedges'] = int(_get_tag(beam,0x300A,0x00D0,0))
        #(300A,00E0) = Number of Compensators [Required (1), integer string]
        bmdata['compensators'] = int(_get_tag(beam,0x300A,0x00E0,0))
        #(300A,00ED) = Number of Boli [Required (1), integer string]
        bmdata['bolus'] = int(_get_tag(beam,0x300A,0x00ED,0))
        #(300A,00F0) = Number of Blocks [Required (1), integer string]
        bmdata['blocks'] = int(_get_tag(beam,0x300A,0x00F0,0))
        #(300A,0304) = Radiation Atomic Number [Conditional (1C), integer string]
        bmdata['atomic_number'] = _get_tag(beam,0x300A,0x0304)
        #(300A,0308) = Scan Mode [Required (1), code string]
        bmdata['scan_mode'] = _get_tag(beam,0x300A,0x0308)
        if bmdata['scan_mode'] == 'MODULATED_SPEC':
            #(300A,0309) = Modulated Scan Mode Type [Conditional (1C), code string]
            bmdata['scan_mode'] = _get_tag(beam,0x300A,0x0309)
        #(300A,0312) = Number of Range Shifters [Required (1), integer string]
        bmdata['range_shifters'] = _get_tag(beam,0x300A,0x0312,0)
        #(300A,0330) = Number of Lateral Spreading Devices [Required (1), integer string]
        bmdata['lateral_spreading'] = _get_tag(beam,0x300A,0x0330,0)
        #(300A,0340) = Number of Range Modulators [Required (1), integer string]
        bmdata['range_modulators'] = _get_tag(beam,0x300A,0x0340,0)
        #(300A,0107) = Applicator Sequence [Optional (3), sequence]
        #(300A,0109) = Applicator Type [Required (1), code string]
        bmdata['applicator'] = _get_tag(_get_tag(beam,0x300A,0x0107,[0])[0],0x300A,0x0109)
        #(300A,0110) = Number of Control Points [Required (1), integer string]
        bmdata['control_points'] = int(_get_tag(beam,0x300A,0x0110,0))
        initial = True
        bmdata['table_lat_dynamic'] = 0
        bmdata['table_long_dynamic'] = 0
        bmdata['table_vert_dynamic'] = 0
        bmdata['table_pitch_dynamic'] = 0
        bmdata['table_roll_dynamic'] = 0
        bmdata['table_rotation_dynamic'] = 0
        bmdata['eccentric_angle_dynamic'] = 0
        bmdata['beam_energy_dynamic'] = 0
        bmdata['dose_rate_dynamic'] = 0
        bmdata['jaw_dynamic'] = None
        bmdata['MLC_dynamic'] = None
        arc_degrees = 0
        #(300A,0111) = Control Point Sequence [Required (1), sequence]
        #(300A,03A8) = Ion Control Point Sequence [Required (1), sequence]
        for cntrlpts in list(_get_tag(beam,0x300A,0x0111,[]))+list(_get_tag(beam,0x300A,0x03A8,[])):
            if initial:
                #(300A,0114) = Nominal Beam Energy [Optional (3), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0114):
                    bmdata['beam_energy'] = float(_get_tag(cntrlpts,0x300A,0x0114))
                else:
                    bmdata['beam_energy'] = None
                #(300A,0115) = Dose Rate Set [Optional (3), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0115):
                    bmdata['dose_rate'] = float(_get_tag(cntrlpts,0x300A,0x0115))
                else:
                    bmdata['dose_rate'] = None
                #(300A,011A) = Beam Limiting Device Position Sequence [Required (1), sequence]
                jaws = dict.fromkeys(['X','Y','ASYMX','ASYMY','MLCX','MLCY'], None)
                for beam_limit in _get_tag(cntrlpts,0x300A,0x011A,[]):
                    #(300A,00B8) = RT Beam Limiting Device [Required (1), code string]
                    limiting_device = _get_tag(beam_limit,0x300A,0x00B8)
                    #(300A,011C) = Leaf/Jaw Positions [Required (1), decimal string]
                    positions = _get_tag(beam_limit,0x300A,0x011C)
                    jaws[limiting_device] = positions
                    if limiting_device in ['MLCX','MLCY']:
                        bmdata['MLC_dynamic'] = 0
                    else:
                        bmdata['jaw_dynamic'] = 0
                #(300A,011E) = Gantry Angle [Conditional (1C), decimal string]
                gantry_angle = float(_get_tag(cntrlpts,0x300A,0x011E,0))
                bmdata['gantry_angle'] = gantry_angle
                #(300A,011F) = Gantry Rotation Direction [Conditional (1C), code string]
                bmdata['gantry_direction'] = _get_tag(cntrlpts,0x300A,0x011F,'NONE')
                #(300A,012A) = Table Top Lateral Position [Conditional (2C), decimal string (mm)]
                if _get_tag(cntrlpts,0x300A,0x012A):
                    bmdata['table_lat'] = float(_get_tag(cntrlpts,0x300A,0x012A))
                else:
                    bmdata['table_lat'] = 0.0
                #(300A,0129) = Table Top Longitudinal Position [Conditional (2C), decimal string (mm)]
                if _get_tag(cntrlpts,0x300A,0x0129):
                    bmdata['table_long'] = float(_get_tag(cntrlpts,0x300A,0x0129))
                else:
                    bmdata['table_long'] = 0.0
                #(300A,0128) = Table Top Vertical Position [Conditional (2C), decimal string (mm)]
                if _get_tag(cntrlpts,0x300A,0x0128):
                    bmdata['table_vert'] = float(_get_tag(cntrlpts,0x300A,0x0128))
                else:
                    bmdata['table_vert'] = 0.0
                #(300A,0140) = Table Top Pitch Angle [Conditional (1C), single (degrees)]
                bmdata['table_pitch'] = float(_get_tag(cntrlpts,0x300A,0x0140,0))
                #(300A,0142) = Table Top Pitch Rotation Direction [Conditional (1C), code string]
                if _get_tag(cntrlpts,0x300A,0x0142,'') == 'CW':
                    bmdata['table_pitch'] *= -1
                #(300A,0144) = Table Top Roll Angle [Conditional (1C), single (degrees)]
                bmdata['table_roll'] = float(_get_tag(cntrlpts,0x300A,0x0144,0)) # table_yaw
                #(300A,0146) = Table Top Roll Rotation Direction [Conditional (1C), code string]
                if _get_tag(cntrlpts,0x300A,0x0146,'') == 'CW':
                    bmdata['table_roll'] *= -1
                #(300A,0122) = Patient Support Angle [Conditional (1C), decimal string]
                bmdata['table_rotation'] = float(_get_tag(cntrlpts,0x300A,0x0122,0))
                #(300A,0123) = Patient Support Rotation Direction [Conditional (1C), code string]
                if _get_tag(cntrlpts,0x300A,0x0123,'') == 'CW':
                    bmdata['table_rotation'] *= -1
                #(300A,0125) = Table Top Eccentric Angle [Conditional (1C), decimal string]
                bmdata['eccentric_angle'] = float(_get_tag(cntrlpts,0x300A,0x0125,0))
                #(300A,0126) = Table Top Eccentric Rotation Direction [Conditional (1C), code string]
                if _get_tag(cntrlpts,0x300A,0x0126,'') == 'CW':
                    bmdata['eccentric_angle'] *= -1
                initial = False
            else:
                #(300A,011A) = Beam Limiting Device Position Sequence [Conditional (1), sequence]
                for beam_limit in _get_tag(cntrlpts,0x300A,0x011A,[]):
                    #(300A,00B8) = RT Beam Limiting Device [Required (1), code string]
                    limiting_device = _get_tag(beam_limit,0x300A,0x00B8)
                    #(300A,011C) = Leaf/Jaw Positions [Required (1), decimal string]
                    positions = _get_tag(beam_limit,0x300A,0x011C)
                    if jaws[limiting_device] != positions:
                        if limiting_device in ['MLCX','MLCY']:
                            bmdata['MLC_dynamic'] = 1
                        elif limiting_device in ['X','Y','ASYMX','ASYMY']:
                            bmdata['jaw_dynamic'] = 1
                #(300A,011E) = Gantry Angle [Conditional (1C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x011E) is not None:
                    gantry_angle_new = float(_get_tag(cntrlpts,0x300A,0x011E))
                    angle_diff = abs(gantry_angle - gantry_angle_new)
                    if angle_diff > 300:
                        angle_diff = 360 - angle_diff
                    arc_degrees += angle_diff
                    gantry_angle = gantry_angle_new
                    #(300A,011F) = Gantry Rotation Direction [Conditional (1C), code string] 
                    gantry_dir = _get_tag(cntrlpts,0x300A,0x011F,'')
                    if gantry_dir != 'NONE' and bmdata['gantry_direction'] != 'NONE' and gantry_dir != bmdata['gantry_direction']:
                        bmdata['gantry_direction'] = 'MIXED'
                #(300A,0128) = Table Top Vertical Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0128) is not None:
                    if bmdata['table_vert'] != float(_get_tag(cntrlpts,0x300A,0x0128)):
                        bmdata['table_vert_dynamic'] = 1
                #(300A,0129) = Table Top Longitudinal Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0129) is not None:
                    if bmdata['table_long'] != float(_get_tag(cntrlpts,0x300A,0x0129)):
                        bmdata['table_long_dynamic'] = 1
                #(300A,012A) = Table Top Lateral Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x012A) is not None:
                    if bmdata['table_lat'] != float(_get_tag(cntrlpts,0x300A,0x012A)):
                        bmdata['table_lat_dynamic'] = 1
                #(300A,0140) = Table Top Pitch Angle [Conditional (1C), single]
                if _get_tag(cntrlpts,0x300A,0x0140) is not None:
                    #(300A,0142) = Table Top Pitch Rotation Direction [Conditional (1C), code string]
                    direction = 1 - 2*(_get_tag(cntrlpts,0x300A,0x0142,'') == 'CW')
                    if bmdata['table_pitch'] != direction * float(_get_tag(cntrlpts,0x300A,0x0140)):
                        bmdata['table_pitch_dynamic'] = 1
                #(300A,0144) = Table Top Roll Angle [Conditional (1C), single]
                if _get_tag(cntrlpts,0x300A,0x0144) is not None:
                    #(300A,0146) = Table Top Roll Rotation Direction [Conditional (1C), code string]
                    direction = 1 - 2*(_get_tag(cntrlpts,0x300A,0x0146,'') == 'CW')
                    if bmdata['table_roll'] != direction * float(_get_tag(cntrlpts,0x300A,0x0144)):
                        bmdata['table_roll_dynamic'] = 1
                #(300A,0122) = Patient Support Angle [Conditional (1C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0122) is not None:
                    #(300A,0123) = Patient Support Rotation Direction [Conditional (1C), code string]
                    direction = 1 - 2*(_get_tag(cntrlpts,0x300A,0x0123,'') == 'CW')
                    if bmdata['table_rotation'] != direction * float(_get_tag(cntrlpts,0x300A,0x0122)):
                        bmdata['table_rotation_dynamic'] = 1
                #(300A,0125) = Table Top Eccentric Angle [Conditional (1C), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0125) is not None:
                    #(300A,0126) = Patient Support Rotation Direction [Conditional (1C), code string]
                    direction = 1 - 2*(_get_tag(cntrlpts,0x300A,0x0126,'') == 'CW')
                    if bmdata['eccentric_angle'] != direction * float(_get_tag(cntrlpts,0x300A,0x0125)):
                        bmdata['eccentric_angle_dynamic'] = 1
                #(300A,0114) = Nominal Beam Energy [Optional (3), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0114):
                    if bmdata['beam_energy'] != float(_get_tag(cntrlpts,0x300A,0x0114)):
                        bmdata['beam_energy_dynamic'] = 1
                #(300A,0115) = Dose Rate Set [Optional (3), decimal string]
                if _get_tag(cntrlpts,0x300A,0x0115):
                    if bmdata['dose_rate'] != float(_get_tag(cntrlpts,0x300A,0x0115)):
                        bmdata['dose_rate_dynamic'] = 1
        bmdata['arc_degrees'] = arc_degrees
        dcm['beams'][beam_ID] = bmdata
    #------- RT Plan : RT Fraction Scheme Module ---------------------------------------------------
    # extract data from fraction group(s) [one or more per treatment plan]
    total_fxs_rx = 0
    brachy_fxs = []
    dcm['plan']['dose'] = None
    dcm['plan']['dose_rx'] = None
    #(300A,0070) = Fraction Group Sequence [Required (1), sequence]
    for fx in _get_tag(plan_data,0x300A,0x0070,[]):
        fxdata = {'planID': plan_ID, 'beams': [], 'modality': None}
        fx_dose = 0
        #(300A,0071) = Fraction Group Number [Required (1), integer string]
        fx_ID = _get_tag(fx,0x300A,0x0071)
        fxdata['fx_groupID'] = fx_ID
        #(300A,0072) = Fraction Group Description [Optional (3), long string]
        fxdata['fx_description'] = _get_tag(fx,0x300A,0x0072)
        #(300A,0078) = Number of Fractions Planned [Required (2), integer string]
        fxdata['num_fxs_rx'] = int(float(_get_tag(fx,0x300A,0x0078,0)))
        total_fxs_rx += fxdata['num_fxs_rx']
        #(300A,0079) = Number of Fraction Pattern Digits per Day [Optional (3), integer string]
        fxdata['pattern_digits_per_day'] = _get_tag(fx,0x300A,0x0079)
        #(300A,007A) = Repeat Fraction Cycle Length [Optional (3), integer string]
        fxdata['repeat_fx_cycle_length'] = _get_tag(fx,0x300A,0x007A)
        #(300A,007B) = Fraction Pattern [Optional (3), long text]
        fxdata['fx_pattern'] = _get_tag(fx,0x300A,0x007B)
        #(300A,0080) = Number of Beams [Required (1), integer string]
        fxdata['num_beams'] = int(_get_tag(fx,0x300A,0x0080))
        #(300A,00A0) = Number of Brachy Application Setups [Required (1), integer string]
        fxdata['num_brachy'] = int(_get_tag(fx,0x300A,0x00A0))
        #(300C,0004) = Referenced Beam Sequence [Conditional (1C), sequence]
        for beam in _get_tag(fx,0x300C,0x0004,[]):
            #(300C,0006) = Referenced Beam Number [Required (1), integer string]
            beam_ID = _get_tag(beam,0x300C,0x0006)
            fxdata['beams'].append(beam_ID)     
            if beam_ID not in dcm['beams']:
                logger.warning(f'Beam ({beam_ID}) improperly referenced in plan ({plan_ID})')
                continue
            elif dcm['beams'][beam_ID]['delivery_type'] not in ['TREATMENT','CONTINUATION']:
                continue
            if not fxdata['modality']:
                fxdata['modality'] = dcm['beams'][beam_ID]['modality']
            elif fxdata['modality'] != dcm['beams'][beam_ID]['modality']:
                fxdata['modality'] = 'MIXED'
            #(300A,008B) = Beam Dose Meaning
            #(300A,0086) = Beam Meterset [Optional (3), decimal string]
            if not dcm['beams'][beam_ID]['cumulative_meterset'] and _get_tag(beam,0x300A,0x0086):
                dcm['beams'][beam_ID]['cumulative_meterset'] = float(_get_tag(beam,0x300A,0x0086))
            #(300A,0090) = Beam Dose Type [Conditional (1C), code string]
            beam_dose_type1 = _get_tag(beam,0x300A,0x0090)
            #(300A,0084) = Beam Dose [Optional (3), decimal string]
            beam_dose1 = float(_get_tag(beam,0x300A,0x0084,0))*100
            #(300A,0092) = Alternative Beam Dose Type [Conditional (1C), code string]
            beam_dose_type2 = _get_tag(beam,0x300A,0x0092)
            #(300A,0091) = Alternative Beam Dose [Optional (3), decimal string]
            beam_dose2 = float(_get_tag(beam,0x300A,0x0091,0))*100
            if beam_dose_type1 == 'PHYSICAL' or beam_dose_type2 == 'EFFECTIVE':
                physical_dose = beam_dose1
                beam_dose = beam_dose1
                effective_dose = beam_dose2
            elif beam_dose_type1 == 'EFFECTIVE' or beam_dose_type2 == 'PHYSICAL':
                physical_dose = beam_dose2
                beam_dose = beam_dose2
                effective_dose = beam_dose1
            else:
                physical_dose = None
                beam_dose = beam_dose1 or beam_dose2
                effective_dose = None
            if beam_dose:
                fx_dose += beam_dose*fxdata['num_fxs_rx']
            dcm['beams'][beam_ID]['fractions'][fx_ID] = {
                    'fx_groupID': fx_ID,
                    'phys_dose_cGy': physical_dose,
                    'eff_dose_cGy': effective_dose
                }
        #(300C,000A) = Referenced Brachy Application Setup Sequence [Conditional (1C), sequence]
        brachy_dose_tot = 0
        for brachy in _get_tag(fx,0x300C,0x000A,[]):
            #(300C,000C) = Referenced Brachy Application Setup Number [Required (1), integer string]
            brachy_ID = _get_tag(brachy,0x300C,0x000C)
            #(300C,00A4) = Brachy Application Setup Dose [Optional (3), decimal string (Gy)]
            brachy_dose = float(_get_tag(brachy,0x300A,0x00A4,0))*100
            if brachy_dose and brachy_dose != brachy_dose_tot:
                brachy_dose_tot += brachy_dose
            #!!! not handling cases where multiple fx groups have bunch of different brachy details...
            brachy_fxs.append({'fx_groupID': fx_ID,
                                    'brachyID': brachy_ID,
                                    'brachy_dose_cGy': brachy_dose})
            if not fxdata['modality']:
                fxdata['modality'] = 'BRACHY'
            elif fxdata['modality'] != 'BRACHY':
                fxdata['modality'] = 'MIXED'
        fx_dose += brachy_dose_tot*fxdata['num_fxs_rx']
        fxdata['rx'] = {}
        fxdata['fx_dose_rx'] = None
        rxdose_max = 0
        #(300C,0050) = Referenced Dose Reference Sequence [Optional (3), sequence]
        for rx in _get_tag(fx,0x300C,0x0050,[]):
            #(300C,0051) = Referenced Dose Reference Number [Required (1), integer string]
            # uniquely identifies Dose Reference specified by Dose Reference Number (300A,0012) in Dose Reference Sequence (300A,0010) in RT Prescription Module.
            rx_ID = _get_tag(rx,0x300C,0x0051)
            if rx_ID not in fxdata['rx']:
                fxdata['rx'][rx_ID] = {'rxID': rx_ID}
            #(300A,0026) = Target Prescription Dose [Optional (3), decimal string (Gy)]
            target_dose_rx = float(_get_tag(rx,0x300A,0x0026,0))*100
            fxdata['rx'][rx_ID]['target_dose_rx'] = target_dose_rx
            if target_dose_rx > rxdose_max:
                rxdose_max = target_dose_rx
        fxdata['fx_dose'] = fx_dose
        dcm['fxs'][fx_ID] = fxdata
        # TO DO???  # can get corresponding ROI_num (3006,0084) from (300A,0010) where (300A,0012) matches rx_id
        if fx_dose and not dcm['plan']['dose']:
            dcm['plan']['dose'] = fx_dose
        elif fx_dose:
            dcm['plan']['dose'] += fx_dose
        if rxdose_max:
            fxdata['fx_dose_rx'] = rxdose_max
            if not dcm['plan']['dose_rx']:
                dcm['plan']['dose_rx'] = rxdose_max
            else:
                dcm['plan']['dose_rx'] += rxdose_max
    dcm['plan']['fxs_rx'] = total_fxs_rx    
    #------- RT Plan : RT Patient Setup Module -----------------------------------------------------
    # extract patient fixation, shielding, positioning, motion-compensation, and setup data                 
    fixation_types = {'CAST': 'casts', 'MOLD': 'molds', 'BITE_BLOCK': 'biteblocks',
                        'RECTAL_BALLOON': 'rectalballoons', 'VACUUM_MOLD': 'vaclocks',
                        'BODY_FRAME': 'bodyframes', 'WHOLE_BODY_POD': 'bodypods', 
                        'HEADREST': 'headrests', 'HEAD_FRAME': 'headframes', 
                        'BREAST_BOARD': 'breastboards', 'MASK': 'masks'}
    shield_types = {'GONAD': 'gonadshields', 'EYE': 'eyeshields', 'GUM': 'gumshields'}
    #(300A,0180) = Patient Setup Sequence [Required (1), sequence]
    posdata = []
    for setup in _get_tag(plan_data,0x300A,0x0180,[]):
        #(300A,0182) = Patient Setup Number [Required (1), integer string]
        setupdata = {'setupID': _get_tag(setup,0x300A,0x0182)}
        #(0018,5100) = Patient Position [Conditional (1C), code string]
        #(300A,0184) = Patient Additional Position [Conditional (1C), long string]
        setupdata['patient_position'] = _get_tag(setup,0x0018,0x5100) or _get_tag(setup,0x300A,0x0184)
        #(300A,01B0) = Setup Technique [Optional (3), code string]
        setupdata['setup_technique'] = _get_tag(setup,0x300A,0x01B0,'')
        #(300A,01B2) = Setup Technique Description [Optional (3), short text]
        setupdata['setup_technique_description'] = _get_tag(setup,0x300A,0x01B2,'')
        for fix_type in fixation_types.values():
            exec(fix_type + ' = []')
        #(300A,0190) = Fixation Device Sequence [Optional (3), sequence]
        for fix in _get_tag(setup,0x300A,0x0190,[]):
            #(300A,0192) = Fixation Device Type [Required (1), code string]
            fix_type = _get_tag(fix,0x300A,0x0192)
            #(300A,0194) = Fixation Device Label [Required (2), short string]
            #(300A,0196) = Fixation Device Description [Optional (3), short text]
            exec(fixation_types[fix_type] + ".append('" + 
                    _get_tag(fix,0x300A,0x0194) + _get_tag(fix,0x300A,0x0196) + "')")
        for fix_type in fixation_types.values():
            exec("setupdata['" + fix_type + "'] = " + str(eval(fix_type)))
        for shield_type in shield_types.values():
            exec(shield_type + ' = []')
        #(300A,01A0) = Shielding Device Sequence [Optional (3), sequence]
        for shield in _get_tag(setup,0x300A,0x01A0,[]):
            #(300A,01A2) = Shielding Device Type [Required (1), code string]
            shield_type = _get_tag(shield,0x300A,0x01A2)
            #(300A,01A4) = Shielding Device Label [Required (2), short string]
            #(300A,01A6) = Shielding Device Descriptioin [Optional (3), short text]
            exec(shield_types[shield_type] + ".append('" + _get_tag(shield,0x300A,0x01A4)
                                + _get_tag(shield,0x300A,0x01A6) + "')")
        for shield_type in shield_types.values():
            exec("setupdata['" + shield_type + "'] = " + str(eval(shield_type)))
        setupdata['motion_compensation_technique'] = []
        setupdata['respiratory_signal'] = []
        #(300A,0410) = Motion Synchronization Sequence [Optional (3), sequence]
        for motion in _get_tag(setup,0x300A,0x0410,[]):
            #(0018,9170) = Respiratory Motion Compensation Technique [Required (1), code string]
            setupdata['motion_compensation_technique'].append(_get_tag(motion,0x0018,0x9170))
            #(0018,9171) = Respiratory Signal Source [Required (1), code string]
            setupdata['respiratory_signal'].append(_get_tag(motion,0x0018,0x9171))        
        posdata.append(setupdata)
    dcm['plan']['setup'] = posdata
    #------- RT Plan : RT Brachy Application Setup Module ------------------------------------------
    # Extract brachytherapy plan data                   
    if len(brachy_fxs) and not dcm['plan']['modality']:
        dcm['plan']['modality'] = 'BRACHY'
    elif len(brachy_fxs):
        dcm['plan']['modality'] = 'MIXED'
    #(300A,0200) = Brachy Treatment Technique [Required (1), code string]
    brachy_technique = _get_tag(plan_data,0x300A,0x0200)
    #(300A,0202) = Brachy Treatment Type [Required (1), code string]
    brachy_type = _get_tag(plan_data,0x300A,0x0202)
    source_data = {}    
    #(300A,0210) = Source Sequence [Required (1), sequence]
    common_isotopes = ['I-125','Ir-192','Pd-103','Cs-131','Cs-137','Co-60','Ru-106','Ra-226']
    for source in _get_tag(plan_data,0x300A,0x0210,[]):
        #(300A,0212) = Source Number [Required (1), integer string]
        source_num = _get_tag(source,0x300A,0x0212)
        srcdata = {}
        #(300A,0214) = Source Type [Required (1), code string]
        srcdata['source_type'] = _get_tag(source,0x300A,0x0214)
        #(300A,0226) = Source Isotope Name [Required (1), long string]
        srcdata['isotope'] = _get_tag(source,0x300A,0x0226)
        # Standardize the source isotope name by finding a match among common isotopes
        srcdata['isotope'] = next((i for i in common_isotopes if i in srcdata['isotope']),'OTHER')
        if srcdata['isotope'] == 'OTHER':
            logger.warning(f"Unexpected source isotope value: {srcdata['isotope']}")
        #(300A,0228) = Source Isotope Half Life [Required (1), decimal string (days)]
        srcdata['halflife'] = float(_get_tag(source,0x300A,0x0228,0))
        #(300A,022A) = Reference Air Kerma Rate [Required (1), decimal string (uGy/h @ 1m)]
        srcdata['air_kerma_rate'] = float(_get_tag(source,0x300A,0x022A,0))
        #(300A,022B) = Source Strength [Conditional (1C), decimal string]
        srcdata['source_strength'] = float(_get_tag(source,0x300A,0x022B,0))
        #(300A,0229) = Source Strength Units [Conditional (1C), code string]
        srcdata['strength_units'] = _get_tag(source,0x300A,0x0229)
        #(300A,022C) = Source Strength Reference Date [Required (1), date (YYYYMMDD)]
        date = _get_tag(source,0x300A,0x022C)
        srcdata['strength_ref_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
        #(300A,022E) = Source Strength Reference Time [Required (1), time (HHMMSS.FFFFFF&ZZXX)]
        time = (_get_tag(source,0x300A,0x022E,'')[:6]+'000000')[:6]
        srcdata['strength_ref_time'] = time[:2]+':'+time[2:4]+':'+time[4:6]
        source_data[source_num] = srcdata
    #(300A,0230) = Application Setup Sequence [Required (1), sequence]
    setups = _get_tag(plan_data,0x300A,0x0230,[])
    device_types = {'MOLD': 'accessory_molds', 'PLAQUE': 'accessory_plaques',  
                    'SHIELD': 'accessory_shields', 'DILATATION': 'accessory_dilatations', 
                    'FLAB': 'accessory_flabs'}
    for fxb in brachy_fxs:
        fxID = fxb['fx_groupID']
        if fxID not in dcm['brachy']:
            dcm['brachy'][fxID] = {'planID': plan_ID, 'fx_groupID': fxID, 
                                    'brachy_technique': brachy_technique,
                                    'brachy_type': brachy_type,
                                    'brachy_dose_cGy': fxb['brachy_dose_cGy'],
                                    'isotope': None, 'source_type': None, 'source_strength': 0,
                                    'strength_units': None, 'air_kerma_rate': 0,
                                    'halflife': 0, 'num_sources': 0, 'num_channels': 0,
                                    'num_control_points': 0, 'total_time': 0,
                                    'total_kerma': 0, 'total_applicator_length': 0,
                                    'application_setup_type': None, 'source_movement': None,
                                    'source_applicator_type': None, 'step_size': None,
                                    'pulse_interval': None, 'accessory_shields': 0,
                                    'accessory_dilatations': 0, 'accessory_molds': 0,
                                    'accessory_plaques': 0, 'accessory_flabs': 0,
                                    'num_application_setups': 0}
        fx_sources = {}
        for setup in setups:
            #(300A,0234) = Application Setup Number [Required (1), integer string]
            setup_num = _get_tag(setup,0x300A,0x0234)
            if setup_num != fxb['brachyID']:
                continue
            dcm['brachy'][fxID]['num_application_setups'] += 1
            #(300A,0232) = Application Setup Type [Required (1), code string]
            app_setup_type = _get_tag(setup,0x300A,0x0232)
            if dcm['brachy'][fxID]['application_setup_type'] is None:
                dcm['brachy'][fxID]['application_setup_type'] = app_setup_type
            elif app_setup_type and app_setup_type != dcm['brachy'][fxID]['application_setup_type']:
                dcm['brachy'][fxID]['application_setup_type'] = 'MULTIPLE'
            #(300A,0242) = Template Type [Optional (3), short string]
            dcm['brachy'][fxID]['template_type'] = _get_tag(setup,0x300A,0x0242)
            #(300A,0250) = Total Reference Air Kerma [Required (1), decimal string]
            dcm['brachy'][fxID]['total_kerma'] += float(_get_tag(setup,0x300A,0x0250,0))
            #(300A,0280) = Channel Sequence [Required (1), sequence]
            channels = _get_tag(setup,0x300A,0x0280,[])
            dcm['brachy'][fxID]['num_channels'] += len(channels)
            for channel in channels:
                #(300A,0110) = Number of Control Points [Required (1), integer string]
                dcm['brachy'][fxID]['num_control_points'] += int(_get_tag(channel,0x300A,0x0110,0))
                #(300A,0286) = Channel Total Time [Required (1), decimal string (sec)]
                dcm['brachy'][fxID]['total_time'] += float(_get_tag(channel,0x300A,0x0286,0))
                #(300A,0288) = Source Movement Type [Required (1), code string]
                src_movement = _get_tag(channel,0x300A,0x0288)
                if dcm['brachy'][fxID]['source_movement'] is None:
                    dcm['brachy'][fxID]['source_movement'] = src_movement
                elif src_movement and src_movement != dcm['brachy'][fxID]['source_movement']:
                    dcm['brachy'][fxID]['source_movement'] = 'MULTIPLE'
                #(300A,028A) = Number of Pulses [Conditional (1C), integer string]
                #(300A,028C) = Pulse Repetition Interval [Conditional (1C), decimal string (sec)]
                if _get_tag(channel,0x300A,0x028C):
                    dcm['brachy'][fxID]['pulse_interval'] = float(_get_tag(channel,0x300A,0x028C))
                #(300A,0296) = Source Applicator Length [Conditional (1C), decimal string (mm)]
                dcm['brachy'][fxID]['total_applicator_length'] += float(_get_tag(channel,0x300A,0x0296,0))
                #(300A,0290) = Source Applicator Number [Optional (3), integer string]
                #(300A,0294) = Source Applicator Name [Optional (3), long string]
                #(300A,0292) = Source Applicator Type [Conditional (1C), code string]
                src_app_type = _get_tag(channel,0x300A,0x0292)
                if dcm['brachy'][fxID]['source_applicator_type'] is None:
                    dcm['brachy'][fxID]['source_applicator_type'] = src_app_type
                elif src_app_type and src_app_type != dcm['brachy'][fxID]['source_applicator_type']:
                    dcm['brachy'][fxID]['source_applicator_type'] = 'MIXED'
                #(300C,000E) = Referenced Source Number [Required (1), integer string]
                source_ref = _get_tag(channel,0x300C,0x000E)
                fx_sources[source_ref] = {
                                    'strength': source_data[source_ref]['source_strength'],
                                    'units': source_data[source_ref]['strength_units']
                                }
                if dcm['brachy'][fxID]['source_type'] is None:
                    dcm['brachy'][fxID]['source_type'] = source_data[source_ref]['source_type']
                elif source_data[source_ref]['source_type'] != dcm['brachy'][fxID]['source_type']:
                    dcm['brachy'][fxID]['source_type'] = 'MULTIPLE'
                if dcm['brachy'][fxID]['isotope'] is None:
                    dcm['brachy'][fxID]['isotope'] = source_data[source_ref]['isotope']
                elif source_data[source_ref]['isotope'] != source_data[source_ref]['isotope']:
                    dcm['brachy'][fxID]['isotope'] = 'MULTIPLE'
                if dcm['brachy'][fxID]['halflife'] == 0:
                    dcm['brachy'][fxID]['halflife'] = source_data[source_ref]['halflife']
                elif source_data[source_ref]['halflife'] != dcm['brachy'][fxID]['halflife']:
                    dcm['brachy'][fxID]['halflife'] = None
                #(300A,02A0) = Source Applicator Step Size [Conditional (1C), decimal string (mm)]
                step_size = _get_tag(channel,0x300A,0x02A0)
                if dcm['brachy'][fxID]['step_size'] is None and step_size:
                    dcm['brachy'][fxID]['step_size'] = _get_tag(channel,0x300A,0x02A0)
                elif step_size and step_size != dcm['brachy'][fxID]['step_size']:
                    dcm['brachy'][fxID]['step_size'] = 'MIXED'
            #(300A,0260) = Brachy Accessory Device Sequence [Optional (3), sequence]
            for device in _get_tag(plan_data,0x300A,0x0260,[]):
                #(300A,0264) = Brachy Accessory Device Type [Required (1), code string]
                dcm['brachy'][fxID][device_types[_get_tag(device,0x300A,0x0264)]] += 1
        dcm['brachy'][fxID]['num_sources'] = len(fx_sources.keys())
        strength_units = list(set([s['units'] for s in fx_sources.values()]))
        if len(strength_units) == 1:
            dcm['brachy'][fxID]['strength_units'] = strength_units[0]
            dcm['brachy'][fxID]['source_strength'] = sum(
                                    [s['strength'] for s in fx_sources.values()]
                                ) / dcm['brachy'][fxID]['num_sources']
        else:
            dcm['brachy'][fxID]['strength_units'] = None
            dcm['brachy'][fxID]['source_strength'] = None
    #------- RT Structure Set : Structure Set Module -----------------------------------------------
    structure_template = dict.fromkeys(['raw_name','description','struct_type','generation_algo',
                                        'interpreted_type','volume_cc','min_dose','max_dose',
                                        'mean_dose','rx_dose','struct_type','min_rx_dose','vol_new',
                                        'raw_name','OAR_full_vol_dose',
                                        'dvh','dvh_new','description','volume_cc','generation_algo',
                                        'dose_ref_type','dose_ref_description','max_rx_dose','label',
                                        'constraint_wt','OAR_limit_dose','OAR_max_dose','mean_dose_new',
                                        'underdose_fx','overdose_fx','min_dose_new','max_dose_new'], None)
    if struct_data:
        #(3006,0020) = Structure Set ROI Sequence [Required (1), sequence]
        for struct in _get_tag(struct_data,0x3006,0x0020,[]):
            #(3006,0022) = ROI Number [Required (1), integer string]
            ROI = _get_tag(struct,0x3006,0x0022)
            if ROI not in dcm['structures']:
                dcm['structures'][ROI] = copy.copy(structure_template)
                dcm['structures'][ROI]['structureID'] = ROI
            #(3006,0026) = ROI Name [Required (2), long string]
            dcm['structures'][ROI]['raw_name'] = _get_tag(struct,0x3006,0x0026)
            #(3006,0028) = ROI Description [Optional (3), short text]
            dcm['structures'][ROI]['description'] = _get_tag(struct,0x3006,0x0028)
            #(3006,0036) = ROI Generation Algorithm [Required (2), code string]
            dcm['structures'][ROI]['generation_algo'] = _get_tag(struct,0x3006,0x0036)
            #(3006,002C) = ROI Volume [Optional (3), decimal string]
            dcm['structures'][ROI]['volume_cc'] = float(_get_tag(struct,0x3006,0x002C,0))    
        #------- RT Structure Set : RT ROI Observations Module -----------------------------------------
        #(3006,0080) = RT ROI Observations Sequence [Required (1), sequence]
        for ROI_obs in _get_tag(struct_data,0x3006,0x0080,[]):
            #(3006,0084) = ROI Number [Required (1), integer string]
            ROI = _get_tag(ROI_obs,0x3006,0x0084)
            if ROI not in dcm['structures']:
                dcm['structures'][ROI] = copy.copy(structure_template)
                dcm['structures'][ROI]['structureID'] = ROI
            #(3006,00A4) = RT ROI Interpreted Type [Required (2), code string] 
            dcm['structures'][ROI]['interpreted_type'] = _get_tag(ROI_obs,0x3006,0x00A4)
            #(3006,0085) = ROI Observation Label [Optional (3), short string] 
            dcm['structures'][ROI]['label'] = _get_tag(ROI_obs,0x3006,0x0085)
    #------- RT Dose : RT Dose Module --------------------------------------------------------------
    for dose_data in [pydicom.dcmread(d) for d in RTdose]:
        #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
        dcm['plan']['dose_ID'] = _get_tag(dose_data,0x0008,0x0018)
        #(3004,000A) = Dose Summation Type [Required (1), code string]
        dose_sum_type = _get_tag(dose_data,0x3004,0x000A)
        dcm['plan']['dose_summation_type'] = dose_sum_type
        if dose_sum_type in ['MULTI_PLAN','FRACTION','BEAM','BRACHY','FRACTION_SESSION','BEAM_SESSION','BRACHY_SESSION','CONTROL_POINT']:
            #(300C,0002) = Referenced RT Plan Sequence [Conditional (1C), sequence]
            for refplan in _get_tag(dose_data,0x300C,0x0002,[]):
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                refID = _get_tag(refplan,0x0008,0x1155)
                if refID != plan_ID:
                    dcm['ref_plans'][refID] = {'planID': plan_ID, 
                                                'referenced_planID': refID,
                                                'refrenced_plan_wt': None}
                    continue
                if dose_sum_type in ['FRACTION','BEAM','BRACHY','FRACTION_SESSION','BEAM_SESSION','BRACHY_SESSION','CONTROL_POINT']:
                    #(300C,0020) = Referenced Fraction Group Sequence [Conditional (1C), sequence]
                    fxgroup = _get_tag(refplan,0x300C,0x0020,[])
                    #(300C,0022) = Referenced Fraction Group Number [Required (1), integer string]
                    fxID = _get_tag(fxgroup,0x300C,0x0022)
                    if dose_sum_type in ['FRACTION','FRACTION_SESSION']:
                        logger.warning(f'RT Dose ({RTdose}) calculated for fraction: {fxID}')
                        # will need to process dose for the fraction and ultimately sum across fractions
                    elif dose_sum_type in ['BEAM','BEAM_SESSION','CONTROL_POINT']:
                        #(300C,0004) = Referenced Beam Sequence [Conditional (1C), sequence]
                        beamseq = _get_tag(refplan,0x300C,0x0004,[])
                        if dose_sum_type in ['BEAM','BEAM_SESSION']:
                            #(300C,0006) = Referenced Beam Number [Required (1), integer string]
                            logger.warning('RT Dose ('+RTdose+') calculated for fraction/beam(s): '
                                +','.join([fxID+'/'+get_tag(b,0x300C,0x0006) for b in beamseq]))
                        if dose_sum_type == 'CONTROL_POINT':
                            #(300C,0006) = Referenced Beam Number [Required (1), integer string]
                            #(300C,00F2) = Referenced Control Point Sequence [Conditional (1C), sequence]
                            #(300C,00F4) = Referenced Start Control Point Index [Required (1), integer string]
                            #(300C,00F6) = Referenced Stop Control Point Index [Required (1), integer string]
                            logger.warning('RT Dose ('+RTdose+') calculated for control point(s) in fraction/beam(s): '
                                +','.join([fxID+'/'+_get_tag(b,0x300C,0x0006) for b in beamseq]))
                    else: # dose_sum_type in ['BRACHY','BRACHY_SESSION']:
                        #(300C,000A) = Referenced Brachy Application Setup Sequence [Conditional (1C), sequence]
                        brachyseq = _get_tag(refplan,0x300C,0x000A,[])
                        #(300C,000C) = Referenced Brachy Application Setup Number [Required (1), integer string]
                        logger.warning('RT Dose ('+RTdose+') calculated for fraction/brachy setup(s): '
                            +','.join([fxID+'/'+_get_tag(b,0x300C,0x000C) for b in brachyseq]))
                else: # dose_sum_type == 'MULTI_PLAN':
                    logger.warning(f'RT Dose ({RTdose}) represents a plan sum across two or more plans')
        elif dose_sum_type == 'RECORD':
            #(3008,0030) = Referenced Treatment Record Sequence [Conditional (1C), sequence]
            for record in _get_tag(dose_data,0x3008,0x0030,[]):
                #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
                record_ID = _get_tag(record,0x0008,0x1155)
                if record_ID not in dcm['treatments']:
                    dcm['treatments'][record_ID] = {'treatmentID': record_ID, 'beams': {}}
                else:
                    continue
                # TO DO:!!! find and process the corresponding RTRECORD data in some way???
                #(300C,0006) Uniquely identifies Beam specified by Referenced Beam Number (300C,0006) in Treatment Session Beam Sequence (3008,0020) of RT Beams Session Record Module within RT Beams Treatment Record referenced in the Referenced Treatment Record Sequence (3008,0030) or in Treatment Session Ion Beam Sequence (3008,0021) of RT Ion Beams Session Record Module within RT Ion Beams Treatment Record referenced in the Referenced Treatment Record Sequence (3008,0030).
            continue
        else: # dose_sum_type == 'PLAN':
            #(3004,0004) = Dose Type [Required (1), code string]
            dcm['plan']['dose_type'] = _get_tag(dose_data,0x3004,0x0004)
            #(3004,0005) = Spatial Transform of Dose [Optional (3), code string]
            dcm['plan']['spatial_transform'] = _get_tag(dose_data,0x3004,0x0005,'NONE')
            #------- RT Dose : RT DVH Module ---------------------------------------------------------------
            # Extract dose volume histogram (DVH) data                  
            #(3004,0042) = DVH Normalization Dose Value [Optional (3), decimal string]
            dcm['plan']['norm_dose'] = _get_tag(dose_data,0x3004,0x0042)
            #(3004,0040) = DVH Normalization Point [Optional (3), decimal string]
            # Coordinates (x, y, z) of common DVH normalization point in the Patient-Based Coordinate System described in Section C.7.6.2.1.1 (mm).
            #(3004,0050) = DVH Sequence [Required (1), sequence]
            for struct in _get_tag(dose_data,0x3004,0x0050,[]):
                #(3004,0060) = DVH Referenced ROI Sequence [Required (1), sequence]
                ROI = _get_tag(struct,0x3004,0x0060,[])
                if len(ROI) > 1:
                    logger.info(f"Skipping DVH import for compound ROI list: [{','.join([_get_tag(r,0x3006,0x0084) for r in ROI])}]")
                    continue
                elif len(ROI) == 0:
                    logger.info('Skipping DVH import (no ROI referenced)')
                    continue
                #(3006,0084) = Referenced ROI Number [Required (1), integer string]
                ROI_num = _get_tag(ROI[0],0x3006,0x0084)
                if ROI_num not in dcm['structures']:
                    dcm['structures'][ROI_num] = copy.copy(structure_template)
                    dcm['structures'][ROI_num]['structureID'] = ROI_num
                logger.info(f"Importing DVH for ROI: {ROI_num} (name: {dcm['structures'][ROI_num]['raw_name']}, type: {dcm['structures'][ROI_num]['interpreted_type']})...")
                #(3004,0002) = DVH Dose Units [Required (1), code string]
                dose_units = _get_tag(struct,0x3004,0x0002)
                # (relative dose DVH) --> convert to absolute dose (cGy)
                if dose_units == 'RELATIVE':
                    if dcm['plan']['norm_dose']:
                        scale = dcm['plan']['norm_dose']/100
                    else:
                        scale = 1
                # (absolute dose DVH) --> convert from Gy to cGy
                else:
                    scale = 100
                #(3004,0070) = DVH Minimum Dose [Optional (3), decimal string]
                min_dose = _get_tag(struct,0x3004,0x0070)
                if min_dose:
                    dcm['structures'][ROI_num]['min_dose'] = float(min_dose)*scale
                #(3004,0072) = DVH Maximum Dose [Optional (3), decimal string]
                max_dose = _get_tag(struct,0x3004,0x0072)
                if max_dose:
                    dcm['structures'][ROI_num]['max_dose'] = float(max_dose)*scale
                #(3004,0074) = DVH Mean Dose [Optional (3), decimal string]
                mean_dose = _get_tag(struct,0x3004,0x0074)
                if mean_dose:
                    dcm['structures'][ROI_num]['mean_dose'] = float(mean_dose)*scale
                #(3004,0004) = Dose Type [Required (1), code string]
                dose_type = _get_tag(struct,0x3004,0x0004)
                #(3004,0052) = DVH Dose Scaling [Required (1), decimal string]
                dose_scaling = float(_get_tag(struct,0x3004,0x0052))
                #(3004,0054) = DVH Volume Units [Required (1), code string]
                volume_units = _get_tag(struct,0x3004,0x0054)
                if volume_units not in ['CM3','PERCENT']:
                    logger.warning("DVH for '"+dcm['structures'][ROI_num]['raw_name']+"' specified in "
                            +'neither absolute nor relative volume units (not currently supported)')
                    dcm['structures'][ROI_num]['dvh'] = None
                    continue
                #(3004,0056) = DVH Number of Bins [Required (1), integer string]
                dosebins = int(_get_tag(struct,0x3004,0x0056))
                #(3004,0058) = DVH Data [Required (1), decimal string]
                dvh = np.array([float(d) for d in _get_tag(struct,0x3004,0x0058,[])])
                if len(dvh) < 3:
                    dcm['structures'][ROI_num]['dvh'] = None
                    continue
                if dosebins*2-1 > len(dvh):
                    dosebins = int(len(dvh)/2)
                DVH_dose = np.cumsum(dvh[::2][:dosebins] * dose_scaling * scale)
                DVH_vol = dvh[1::2][:dosebins]
                #(3004,0001) = DVH Type [Required (1), code string]
                dvh_type = _get_tag(struct,0x3004,0x0001)
                if dvh_type == 'CUMULATIVE':
                    if volume_units == 'CM3':
                        total_vol = np.max(DVH_vol)
                        if not total_vol:
                            dcm['structures'][ROI_num]['dvh'] = None
                            continue
                        if not dcm['structures'][ROI_num]['volume_cc']:
                            dcm['structures'][ROI_num]['volume_cc'] = total_vol
                        DVH_vol = 100 * DVH_vol / total_vol
                    else: # volume_units == 'PERCENT'
                        if not dcm['structures'][ROI_num]['mean_dose']:
                            diff_vol = np.insert(np.diff(DVH_vol), 0, 0)
                            dcm['structures'][ROI_num]['mean_dose'] = np.sum(DVH_dose * diff_vol) / np.sum(diff_vol)
                elif dvh_type == 'DIFFERENTIAL':
                    # TO DO: using exported data from Eclipse for differential DVHs, volume was scaled unreliably
                    # (not sure if this is an error/quirk from Eclipse v. other TPS)
                    # !!!!! need to verify behavior in other exported DICOMs from other TPS types
                    volume_scale = 1 # dvh[0]
                    DVH_vol = [v*volume_scale for v in DVH_vol]
                    total_vol = np.sum(DVH_vol)
                    if not total_vol:
                        dcm['structures'][ROI_num]['dvh'] = None
                        continue
                    # absolute volume (cc) DVH --> convert to relative volume (%) DVH
                    if volume_units == 'CM3':
                        if not dcm['structures'][ROI_num]['volume_cc']:
                            dcm['structures'][ROI_num]['volume_cc'] = total_vol
                        if not dcm["structures"][ROI_num]['mean_dose']:
                            dcm['structures'][ROI_num]['mean_dose'] = np.sum(DVH_dose * DVH_vol) / total_vol
                        DVH_vol = (100 * np.cumsum(DVH_vol[::-1]))[::-1] / total_vol
                    # relative volume (%) DVH
                    else: # volume_units == 'PERCENT'
                        if not dcm['structures'][ROI_num]['mean_dose']:
                            dcm['structures'][ROI_num]['mean_dose'] = np.sum(DVH_dose * DVH_vol) / total_vol
                else: # dvh_type = 'NATURAL' [see https://doi.org/10.1118/1.595815]
                    logger.warning("'NATURAL' type is not supported for DVH for '"
                        +dcm['structures'][ROI_num]['raw_name']+"'")
                    dcm['structures'][ROI_num]['dvh'] = None
                    continue
                dmax = float(dcm['structures'][ROI_num]['max_dose'])
                dmin = float(dcm['structures'][ROI_num]['min_dose'])
                DVH_dose = [100*(d-dmin)/(dmax-dmin) for d in DVH_dose]
                dcm['structures'][ROI_num]['dvh'] = _DVH_interp(DVH_dose,DVH_vol,[i/2 for i in range(1,200)])
    # Calculate a new DVH for each structure and extract structure surface area and updated volume
    for dvh in calculate_dvh(structure=struct_data, dose=dose_data, rois=list(dcm['structures'].keys())):
        dcm['structures'][dvh['ROI']]['vol_new'] = dvh['volume']
        dcm['structures'][dvh['ROI']]['vol_incomplete_new'] = dvh['incomplete_dose_volume']        
        dcm['structures'][dvh['ROI']]['min_dose_new'] = float(dvh['dose_min'])
        dcm['structures'][dvh['ROI']]['mean_dose_new'] = float(dvh['dose_mean'])
        dcm['structures'][dvh['ROI']]['max_dose_new'] = float(dvh['dose_max'])
        dcm['structures'][dvh['ROI']]['integral_dose_new'] = dvh['dose_integral']
        dcm['structures'][dvh['ROI']]['surface_area_new'] = dvh['surface_area']
        dcm['structures'][dvh['ROI']]['dvh_new'] = list(dvh['dvh'])
    #------- RT Plan : RT Prescription Module ------------------------------------------------------
    #(300A,000E) = Prescription Description [Optional (3), short text]
    dcm['plan']['rx_description'] = _get_tag(plan_data,0x300A,0x000E)
    #(300A,0010) = Dose Reference Sequence [Optional (3), sequence]
    rxdose_max = 0
    for rx in _get_tag(plan_data,0x300A,0x0010,[]):
        #(300A,0014) = Dose Reference Structure Type [Required (1), code string]
        structreftype = _get_tag(rx,0x300A,0x0014)        
        #(300A,0020) = Dose Reference Type [Required (1), code string]
        dosereftype = _get_tag(rx,0x300A,0x0020)
        #(300A,0026) = Target Prescription Dose [Optional (3), decimal string (Gy)]
        rx_dose = _get_tag(rx,0x300A,0x0026)
        if rx_dose:
            rx_dose = float(rx_dose)*100
            if rx_dose > rxdose_max and dosereftype == 'TARGET' and structreftype != 'COORDINATES':
                rxdose_max = rx_dose
        #(3006,0084) = Referenced ROI Number [Conditional (1C), integer string]
        ROI_num = _get_tag(rx,0x3006,0x0084)
        if not ROI_num:
            continue
        if ROI_num not in dcm['structures'].keys():
            dcm['structures'][ROI_num] = copy.copy(structure_template)
            dcm['structures'][ROI_num]['structureID'] = ROI_num
        dcm['structures'][ROI_num]['struct_type'] = structreftype
        dcm['structures'][ROI_num]['dose_ref_type'] = dosereftype
        #(300A,0012) = Dose Reference Number [Required (1), integer string]
        dcm['structures'][ROI_num]['dose_ID'] = _get_tag(rx,0x300A,0x0012)
        #(300A,0016) = Dose Reference Description [Optional (3), long string]
        dcm['structures'][ROI_num]['dose_ref_description'] = _get_tag(rx,0x300A,0x0016)
        #(300A,0021) = Constraint Weight [Optional (3), decimal string]
        dcm['structures'][ROI_num]['constraint_wt'] = _get_tag(rx,0x300A,0x0021)
        #(300A,0028) = Target Underdose Volume Fraction [Optional (3), decimal string]
        dcm['structures'][ROI_num]['underdose_fx'] = _get_tag(rx,0x300A,0x0028, 0)
        #(300A,002D) = Organ at Risk Overdose Volume Fraction [Optional (3), decimal string]
        dcm['structures'][ROI_num]['overdose_fx'] = _get_tag(rx,0x300A,0x002D)
        #(300A,0025) = Target Minimum Dose [Optional (3), decimal string (Gy)]
        if _get_tag(rx,0x300A,0x0025):
            dcm['structures'][ROI_num]['min_rx_dose'] = float(_get_tag(rx,0x300A,0x0025))*100
        #(300A,0027) = Target Maximum Dose [Optional (3), decimal string (Gy)]
        if _get_tag(rx,0x300A,0x0027):
            dcm['structures'][ROI_num]['max_rx_dose'] = float(_get_tag(rx,0x300A,0x0027))*100
        #(300A,002A) = Organ at Risk Full-volume Dose [Optional (3), decimal string]
        if _get_tag(rx,0x300A,0x002A):
            dcm['structures'][ROI_num]['OAR_full_vol_dose'] = float(_get_tag(rx,0x300A,0x002A))*100
        #(300A,002B) = Organ at Risk Limit Dose [Optional (3), decimal string]
        if _get_tag(rx,0x300A,0x002B):
            dcm['structures'][ROI_num]['OAR_limit_dose'] = float(_get_tag(rx,0x300A,0x002B))*100
        #(300A,002C) = Organ at Risk Maximum Dose [Optional (3), decimal string]
        if _get_tag(rx,0x300A,0x002C):
            dcm['structures'][ROI_num]['OAR_max_dose'] = float(_get_tag(rx,0x300A,0x002C))*100
        if rx_dose:
            dcm['structures'][ROI_num]['rx_dose'] = rx_dose
    if rxdose_max and not dcm['plan']['dose_rx']:
        dcm['plan']['dose_rx'] = rxdose_max
    #------- RT (Ion) Beams / Brachy Treatment Record : General Study Module --------------------------------------
    for record_data in [pydicom.dcmread(r) for r in RTrecord]:
        #(0008,0018) = SOP Instance UID [Required (1), unique identifier]
        recID = _get_tag(record_data,0x0008,0x0018)
        if recID not in dcm['treatments'].keys():
            dcm['treatments'][recID] = {'treatmentID': recID, 'planID': plan_ID, 'beams': [], 'brachy': [], 'measurements': []}
        else:
            continue
        #(0020,0013) = Instance Number [Required (1), integer string]
        instance_num = _get_tag(record_data,0x0020,0x0013)
        dcm['treatments'][recID]['instance_num'] = instance_num
        #(300A,0078) = Number of Fractions Planned [Required (2), integer string]
        dcm['treatments'][recID]['num_fxs_rx'] = _get_tag(record_data,0x300A,0x0078)
        #(300C,0022) = Referenced Fraction Group Number [Optional (3), integer string]
        dcm['treatments'][recID]['fx_groupID'] = _get_tag(record_data,0x300C,0x0022)
        #(3008,0200) = Current Treatment Status [Required (1), code string]
        dcm['treatments'][recID]['current_status'] = _get_tag(record_data,0x3008,0x0200)
        #(3008,0202) = Treatment Status Comment [Optional (3), short text]
        dcm['treatments'][recID]['comment'] = _get_tag(record_data,0x3008,0x0202)
        #(3008,0054) = First Treatment Date [Required (2), date (YYYYMMDD is DICOM standard)]
        date = _get_tag(record_data,0x3008,0x0054)
        if date:
            dcm['treatments'][recID]['first_tx_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
        else:
            dcm['treatments'][recID]['first_tx_date'] = None
        #(3008,0056) = Most Recent Treatment Date [Required (2), date (YYYYMMDD is DICOM standard)]
        date = _get_tag(record_data,0x3008,0x0056)
        if date:
            dcm['treatments'][recID]['most_recent_tx_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
        else:
            dcm['treatments'][recID]['most_recent_tx_date'] = None
        #(3008,0250) = Treatment Date [Required (2), date (YYYYMMDD is DICOM standard)]
        date = _get_tag(record_data,0x3008,0x0250)
        logger.info(f'Processing RT treatment #{instance_num} (ID: {recID}, date: {date})...')
        if date:
            dcm['treatments'][recID]['treatment_date'] = date[:4]+'-'+date[4:6]+'-'+date[6:8]
        else:
            dcm['treatments'][recID]['treatment_date'] = None
        #(3008,0251) = Treatment Time [Required (2), time (HHMMSS.FFFFFF&ZZXX is DICOM standard)]
        time = (_get_tag(record_data,0x3008,0x0251,'')[:6]+'000000')[:6]
        dcm['treatments'][recID]['treatment_time'] = time[:2]+':'+time[2:4]+':'+time[4:6]
        #(3008,0010) = Measured Dose Reference Sequence [Required (1), sequence]
        for measure in _get_tag(record_data,0x3008,0x0010,[]):
            #(3008,0014) = Measured Dose Type [Required (2), code string]
            measure_type = _get_tag(measure,0x3008,0x0014,'')
            #(3008,0012) = Measured Dose Description [Optional (3), short text]
            measure_description = _get_tag(measure,0x3008,0x0012,'')
            #(3004,0002) = Dose Units [Required (1), code string]
            measure_units = _get_tag(measure,0x3004,0x0002)
            #(3008,0016) = Measured Dose Value [Required (2), decimal string]
            measure_value = _get_tag(measure,0x3008,0x0016)
            # TO DO: convert absolute measured dose into relative dose compared with calculated dose reference sequence (3008,0070) (I'm assuming this is a proper thing to do???)
            #if measure_units != 'RELATIVE':
            #   ref_number = _get_tag(measure,0x3008,0x0064)
            #   if not ref_number:
            #       ref_number = _get_tag(measure,0x300C,0x0051)
            dcm['treatments'][recID]['measurements'].append({'type': measure_type,
                                                            'description': measure_description,
                                                            'value': measure_value,
                                                            'units': measure_units})
        #(300A,0206) = Treatment Machine Sequence [Required (1), sequence (length=1)]
        for machine in _get_tag(record_data,0x300A,0x0206,[]):
            #(300A,0206) = Treatment Machine Name [Required (2), short string]
            dcm['treatments'][recID]['machine_name'] = _get_tag(machine,0x300A,0x00B2)
            #(0008,0070) = Manufacturer [Required (2), long string]
            dcm['treatments'][recID]['manufacturer'] = _get_tag(machine,0x0008,0x0070) 
            #(0008,1090) = Manufacturer's Model Name [Required (2), long string]
            dcm['treatments'][recID]['model_name'] = _get_tag(machine,0x0008,0x1090)
            #(0008,0080) = Institution Name [Required (2), long string]
            dcm['treatments'][recID]['institution_name'] = _get_tag(machine,0x0008,0x0080)
            #(0018,1000) = Device Serial Number [Required (2), long string]
            dcm['treatments'][recID]['serial_no'] = _get_tag(machine,0x0018,0x1000)
        #(3008,0020) = Treatment Session Beam Sequence [Required (1), sequence]
        #(3008,0021) = Treatment Session Ion Beam Sequence [Required (1), sequence]
        for beam in list(_get_tag(record_data,0x3008,0x0020,[]))+list(_get_tag(record_data,0x3008,0x0021,[])):
            # Ignore non-treatment beams (e.g. setup fields, port films, verifications)
            #(300A,00CE) = Treatment Delivery Type [Required (2), code string]
            if _get_tag(beam,0x300A,0x00CE,'UNKNOWN') not in ['TREATMENT','CONTINUATION']:
                continue
            #(300C,0006) = Referenced Beam Number [Optional (3), integer string]
            beam_ID = _get_tag(beam,0x300C,0x0006)
            beam_data = {'beamID': beam_ID, 'shift_lat': None, 'shift_long': None,
                            'shift_vert': None, 'adj_pitch': None, 'adj_roll': None, 'adj_yaw': None}
            #(3008,0022) = Current Fraction Number [Required (2), integer string]
            beam_data['fraction_number'] = _get_tag(beam,0x3008,0x0022)
            #(3008,002A) = Treatment Termination Status [Required (1), code string]
            beam_data['termination_status'] = _get_tag(beam,0x3008,0x002A)
            #(3008,002C) = Treatment Verification Status [Required (2), code string]
            beam_data['verification_status'] = _get_tag(beam,0x3008,0x002C)
            tbl_lat = None
            tbl_long = None
            tbl_vert = None
            tbl_pitch = None
            tbl_roll = None
            tbl_rotation = None
            #(3008,0040) = Control Point Delivery Sequence [Required (1), sequence]
            #(3008,0041) = Ion Control Point Delivery Sequence [Required (1), sequence]
            cntrlpts = list(_get_tag(beam,0x3008,0x0040,[]))+list(_get_tag(beam,0x3008,0x0041,[]))
            #(3008,003B) = Delivered Treatment Time [Optional (3), decimal string (sec)]
            beam_data['delivery_time'] = _get_tag(beam,0x3008,0x003B)
            if not beam_data['delivery_time'] and len(cntrlpts) > 1:
                #(3008,0025) = Treatment Control Point Time [Required (1), time]
                beam_data['delivery_time'] = _get_tag(cntrlpts[-1],0x3008,0x0025)-_get_tag(cntrlpts[0],0x3008,0x0025)
            for cntrlpt in cntrlpts:
                #(300A,0128) = Table Top Vertical Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpt,0x300A,0x0128):
                    tbl_vert = float(_get_tag(cntrlpt,0x300A,0x0128))
                #(300A,0129) = Table Top Longitudinal Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpt,0x300A,0x0129):
                    tbl_long = float(_get_tag(cntrlpt,0x300A,0x0129))
                #(300A,012A) = Table Top Lateral Position [Conditional (2C), decimal string]
                if _get_tag(cntrlpt,0x300A,0x012A):
                    tbl_lat = float(_get_tag(cntrlpt,0x300A,0x012A))
                #(300A,0140) = Table Top Pitch Angle [Conditional (1C), single]
                if _get_tag(cntrlpt,0x300A,0x0140) is not None and not tbl_pitch:
                    tbl_pitch = _get_tag(cntrlpt,0x300A,0x0140,0)
                    #(300A,0142) = Table Top Pitch Rotation Direction [Conditional (1C), code string]
                    tbl_pitch_dir = _get_tag(cntrlpt,0x300A,0x0142,'NONE')
                    if tbl_pitch_dir == 'CW':
                        tbl_pitch *= -1
                    elif tbl_pitch_dir == 'NONE':
                        tbl_pitch = 0
                #(300A,0122) = Patient Support Angle [Conditional (1C), decimal string]
                if _get_tag(cntrlpt,0x300A,0x0122) is not None and not tbl_rotation:
                    tbl_rotation = _get_tag(cntrlpt,0x300A,0x0122,0) #  (Table Top Yaw Angle)
                    #(300A,0123) = Patient Support Rotation Direction [Conditional (1C), code string]
                    tbl_rotation_dir = _get_tag(cntrlpt,0x300A,0x0123,'NONE')
                    if tbl_rotation_dir == 'CW':
                        tbl_rotation *= -1
                    elif tbl_rotation_dir == 'NONE':
                        tbl_rotation = 0
                #(300A,0144) = Table Top Roll Angle [Conditional (1C), single]
                if _get_tag(cntrlpt,0x300A,0x0144) is not None and not tbl_roll:
                    tbl_roll = _get_tag(cntrlpt,0x300A,0x0144)
                    #(300A,0146) = Table Top Roll Rotation Direction [Conditional (1C), code string]
                    tbl_roll_dir = _get_tag(cntrlpt,0x300A,0x0146,'NONE')
                    if tbl_roll_dir == 'CW':
                        tbl_roll *= -1
                    elif tbl_roll_dir == 'NONE':
                        tbl_roll = 0
            if beam_ID in dcm['beams'].keys():
                if dcm['beams'][beam_ID]['table_lat'] and tbl_lat:
                    beam_data['shift_lat'] = tbl_lat - float(dcm['beams'][beam_ID]['table_lat'])
                if dcm['beams'][beam_ID]['table_long'] and tbl_long:
                    beam_data['shift_long'] = tbl_long - float(dcm['beams'][beam_ID]['table_long'])
                if dcm['beams'][beam_ID]['table_vert'] and tbl_vert:
                    beam_data['shift_vert'] = tbl_vert - float(dcm['beams'][beam_ID]['table_vert'])
                if dcm['beams'][beam_ID]['table_pitch'] and tbl_pitch:
                    beam_data['adj_pitch'] = tbl_pitch - float(dcm['beams'][beam_ID]['table_pitch'])
                if dcm['beams'][beam_ID]['table_roll'] and tbl_roll:
                    beam_data['adj_roll'] = tbl_roll - float(dcm['beams'][beam_ID]['table_roll'])
                if dcm['beams'][beam_ID]['table_rotation'] and tbl_rotation:
                    beam_data['adj_yaw'] = tbl_rotation - float(dcm['beams'][beam_ID]['table_rotation'])
            dcm['treatments'][recID]['beams'].append(beam_data)
        #(3008,0110) = Treatment Session Application Setup Sequence [Required (1), sequence]
        for brachy in list(_get_tag(record_data,0x3008,0x0110,[])):
            #(300C,000C) = Referenced Brachy Application Setup Number [Optional (3), integer string]
            brachy_ID = _get_tag(brachy,0x300C,0x0006)
            #(300A,0250) = Total Reference Air Kerma [Required (1), decimal string]
            total_kerma = _get_tag(brachy,0x300A,0x0250,0)
            #(3008,002A) = Treatment Termination Status [Required (1), code string]
            termination_status = _get_tag(brachy,0x3008,0x002A,'UNKNOWN')
            #(300A,00CE) = Treatment Delivery Type [Required (2), code string]
            delivery_type = _get_tag(brachy,0x300A,0x00CE,'UNKNOWN')
            #(3008,002C) = Treatment Verification Status [Required (2), code string]
            verification = _get_tag(brachy,0x3008,0x002C,'UNKNOWN')
            #(3008,0022) = Current Fraction Number [Required (2), integer string]
            fx_num = _get_tag(brachy,0x3008,0x0022)
            #(3008,0130) = Recorded Channel Sequence [Required (1), sequence]
            #(3008,0134) = Delivered Channel Total Time [Required (1), decimal string]
            time = sum([_get_tag(c,0x3008,0x0134,0) for c in _get_tag(brachy,0x3008,0x0130,[])])
            brachy_data = {'brachyID': brachy_ID, 'total_air_kerma': total_kerma,
                                'termination_status': termination_status, 'delivery_type': delivery_type,
                                'verification_status': verification, 'fx': fx_num, 'delivery_time': time}
            dcm['treatments'][recID]['brachy'].append(brachy_data)
    return(json.dumps(dcm))

# STRUCTURE COORDINATES MAY NOT AGREE WITH DOSE COORDINATES (DIRECTION COSINES - verify)
# ^^ where does this come from?  no cosines present in structure set itself --> do I need to pull from referenced image(s) somehow?
# TO CONSIDER: STUDY OPTIMIZATION OF RESOLUTION W/ COMPUTE PERFORMANCE V. OUTCOME / ACCURACY -- obviously for smaller structures, need higher res, but this is already handled to some degree (>50 subdivisions in given axis)
def calculate_dvh(structure, dose, rois=[], res=[0.5,0.5,0.5]):
    """
    Calculates Dose-Volume Histograms (DVHs) for specified regions of interest (ROIs) within a radiation therapy structure, using dose grid data.

    This function processes the given DICOM RT structure and dose objects to compute DVHs. It starts by extracting relevant dose 
    and geometric data from the DICOM files, which includes pixel data, dose grid scaling, and image orientation. The function handles 
    varying pixel resolutions and orientations, ensuring accuracy in calculations. It can process multiple ROIs, handling each one based on 
    the provided structure and dose. The function also manages cases of incomplete dose data, issuing warnings for structures outside the 
    dose grid or with missing dose information. For each ROI, it calculates minimum, maximum, mean doses, and the integral dose, along with DVHs. 
    The function is robust, handling various inconsistencies and edge cases, and logs warnings for potential issues.

    Args:
        structure (pydicom Dataset): The DICOM RT structure dataset, containing structure data.
        dose (pydicom Dataset): The DICOM RT dose dataset, containing dose grid data.
        rois (list of str, optional): A list of Region of Interest (ROI) ids to calculate DVHs for. 
                                      Defaults to an empty list, which means all ROIs are processed.
        res (list of float, optional): A 3-element list specifying the resolution [x, y, z] to which the 
                                       dose grid should be resampled. Defaults to [0.5, 0.5, 0.5].

    Returns:
        list of dict: A list of dictionaries, each representing an ROI. Each dictionary contains detailed 
                      DVH information such as minimum, mean, and maximum doses, integral dose, volume, 
                      surface area, and the DVH itself. Incomplete dose volume and warnings are also included 
                      for ROIs with partial or missing dose information.
    """
    #(7FE0,0010) = Pixel Data [Conditional (1C), OB or OW bit stream] 
    #   (NOTE: .pixel_array is a pydicom-specific pixel data extract)
    #(3004,000E) = Dose Grid Scaling [Conditional (1C), decimal string]
    dose_grid = dose.pixel_array * float(_get_tag(dose,0x3004,0x000E,1))
    # Convert dose_grid numpy array to VTK array
    dose_array = vtk_np.numpy_to_vtk(num_array=dose_grid.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    [n_frame,n_row,n_col] = dose_grid.shape
    #(3004,0042) = DVH Normalization Dose Value [Optional (3), decimal string]
    norm_dose = _get_tag(dose,0x3004,0x0042,1)
    #(3004,0040) = DVH Normalization Point [Optional (3), decimal string]
    # Coordinates (x, y, z) of common DVH normalization point in the Patient-Based Coordinate System described in Section C.7.6.2.1.1 (mm).
    norm_point = _get_tag(dose,0x3004,0x0040)
    #(3004,0004) = Dose Type [Required (1), code string]
    dose_type = _get_tag(dose,0x3004,0x0004)
    #(0020,0032) = Image Position (Patient) [Required (1), decimal string], specifies the x, y, z coords of the upper left hand corner of the image; it is the center of the first voxel transmitted
    [orig_x,orig_y,orig_z] = [float(pos) for pos in _get_tag(dose,0x0020,0x0032)]
    #(0020,0037) = Image Orientation (Patient) [Required (1), decimal string]
    # Xx, Xy, Xz = direction cosines for the row
    # Yx, Yy, Yz = direction cosines for the column
    # If Anatomical Orientation Type (0010,2210) is absent or has a value of BIPED, the x-axis is increasing to the left hand side of the patient; the y-axis is increasing to the posterior side of the patient; the z-axis is increasing toward the head of the patient.
    [Xx,Xy,Xz,Yx,Yy,Yz] = [float(dir_cos) for dir_cos in _get_tag(dose,0x0020,0x0037)]
    Zx, Zy, Zz = np.cross([Xx, Xy, Xz], [Yx, Yy, Yz])
    orientation_matrix = vtk.vtkMatrix4x4()
    elements = [Xx,Xy,Xz,Yx,Yy,Yz,Zx,Zy,Zz]
    for i in range(3):
        for j in range(3):
            orientation_matrix.SetElement(i, j, elements[i + j * 3])
    #(0028,0030) = Pixel Spacing [Required (1), decimal string (mm)]
    # row_spacing corresponds to the y-coordinate spacing
    # col_spacing corresponds to the x-coordinate spacing
    [row_spacing,col_spacing] = [float(spc) for spc in _get_tag(dose,0x0028,0x0030,[1,1])]
    #(3004,000C) = Grid Frame Offset Vector [Conditional (1C), decimal string]
    grid_offset = _get_tag(dose,0x3004,0x000C)
    #(0018,0050) = Slice Thickness [Required (2), decimal string (mm)]
    slice_thickness = _get_tag(dose,0x0018,0x0050)
    if slice_thickness and not grid_offset:
        frame_spacing = float(slice_thickness)
    elif grid_offset and len(set([round(grid_offset[i+1]-grid_offset[i],3) for i in range(len(grid_offset)-1)])) == 1:
        frame_spacing = round(grid_offset[1]-grid_offset[0],3)
    else:
        frame_spacing = 0
        warnings.warn('Unspecified or inconsistent frame spacing.')
        return([])
        # QUESTION: HOW TO HANDLE FRAME SPACING OF 0??? WHEN DOES THIS HAPPEN??
    # Create and populate a vtkImageData with dose grid data
    dose_image = vtk.vtkImageData()
    # Set the origin, spacing, and dimensions; note: the VTK image origin defines the physical location of the corner of the first voxel in the image
    dose_grid_min_extent = [round(orig_x-col_spacing*0.5,10),round(orig_y-row_spacing*0.5,10),round(orig_z-frame_spacing*0.5,10)]
    dose_grid_max_extent = [round(orig_x+col_spacing*(n_col-0.5),10),round(orig_y+row_spacing*(n_row-0.5),10),round(orig_z+frame_spacing*(n_frame-0.5),10)]
    dose_image.SetOrigin(dose_grid_min_extent[0],dose_grid_min_extent[1],dose_grid_min_extent[2])
    dose_image.SetDimensions(n_col, n_row, n_frame)
    dose_image.SetSpacing(col_spacing, row_spacing, frame_spacing)
    # Set the scalar data for dose_image
    dose_image.GetPointData().SetScalars(dose_array)
    # Upsample the dose grid to specified resolution
    resample = vtk.vtkImageResample()
    resample.SetInputData(dose_image)
    resample.SetOutputSpacing(min(res[0],col_spacing), min(res[1],row_spacing), min(res[2],frame_spacing))
    resample.SetInterpolationModeToLinear()
    resample.Update()
    resampled_dose = resample.GetOutput()
    res = list(resampled_dose.GetSpacing())
    orig_x -= res[0]*0.5
    orig_y -= res[1]*0.5
    orig_z -= res[2]*0.5
    resampled_dose.SetOrigin(orig_x,orig_y,orig_z)
    # Obtain the new origin, spacing, and dimensions from the resampled image
    [n_col, n_row, n_frame] = resampled_dose.GetDimensions()
    return_data = []
    struct_names = {_get_tag(s,0x3006,0x0022): _get_tag(s,0x3006,0x0026) for s in _get_tag(structure,0x3006,0x0020,[])}
    for roi in _rois_from_rtstruct(structure,rois=rois,res=res):
        mesh_bounds = roi['mesh'].GetOutput().GetBounds()
        # Avoid DVH calculation for structures completely outside of dose grid extents
        if (mesh_bounds[1] < dose_grid_min_extent[0] or mesh_bounds[3] < dose_grid_min_extent[1] or 
                mesh_bounds[5] < dose_grid_min_extent[2] or mesh_bounds[0] > dose_grid_max_extent[0] or 
                mesh_bounds[2] > dose_grid_max_extent[1] or mesh_bounds[4] > dose_grid_max_extent[2]):
            logger.warning(f"Dose not calculated for structure {roi['ROI']} ('{struct_names[roi['ROI']]}')")
            continue
        structure_stencil = _stencil_from_mesh(roi['mesh'].GetOutputPort(),dose_image,orientation_matrix)
        image_stencil = vtk.vtkImageStencil()
        image_stencil.SetInputData(dose_image)
        image_stencil.SetStencilConnection(structure_stencil.GetOutputPort())
        image_stencil.SetBackgroundValue(-1)
        image_stencil.Update()
        # TO DO: get surface min and max doses (interpolated from rough dose grid) across all points on surface mesh --> update min/max doses accordingly
        structure_doses = vtk_np.vtk_to_numpy(image_stencil.GetOutput().GetPointData().GetScalars())
        structure_doses = structure_doses[structure_doses >= 0]
        if len(structure_doses) == 0: # this is test for error Tianjun was encountering ...
            logger.warning(f"Failed DVH calculation for ROI {roi['ROI']} ('{struct_names[roi['ROI']]}') [structure located outside dose grid]")
            continue
        dose_min, dose_max = np.min(structure_doses), np.max(structure_doses)
        incomplete_dose = (mesh_bounds[0] < orig_x or mesh_bounds[2] < orig_y or mesh_bounds[4] < orig_z or 
            mesh_bounds[1] > orig_x+res[0]*n_col or mesh_bounds[3] > orig_y+res[1]*n_row or 
            mesh_bounds[5] > orig_z+res[2]*n_frame)
        incomplete_orig_dose = (mesh_bounds[0] < dose_grid_min_extent[0] or mesh_bounds[2] < dose_grid_min_extent[1] or 
            mesh_bounds[4] < dose_grid_min_extent[2] or mesh_bounds[1] > dose_grid_max_extent[0] or 
            mesh_bounds[3] > dose_grid_max_extent[1] or mesh_bounds[5] > dose_grid_max_extent[2])
        if roi['type'] != 'EXTERNAL' and not incomplete_dose:
            structure_stencil = _stencil_from_mesh(roi['mesh'].GetOutputPort(), resampled_dose, orientation_matrix)
            image_stencil = vtk.vtkImageStencil()
            image_stencil.SetInputData(resampled_dose)
            image_stencil.SetStencilConnection(structure_stencil.GetOutputPort())
            image_stencil.SetBackgroundValue(-1)
            image_stencil.Update()
            structure_doses = vtk_np.vtk_to_numpy(image_stencil.GetOutput().GetPointData().GetScalars())
            structure_doses = structure_doses[structure_doses >= 0]
            dose_min = min(dose_min, np.min(structure_doses))
            dose_max = max(dose_max, np.max(structure_doses))
            dose_mean = np.mean(structure_doses)
            volume_missing = 0
            integral_dose = np.prod(res)*0.001*structure_doses.sum()
        # Reduce compute cost for body ('EXTERNAL') contour by stenciling over original (not upsampled) dose grid
        else: # ROI_type == 'EXTERNAL' or incomplete_dose
            # Flag DVH calculation for EXTERNAL structures without complete dose grid data
            if incomplete_orig_dose:
                volume_missing = max(0, roi['volume']-len(structure_doses)*col_spacing*row_spacing*frame_spacing*0.001)
                logger.warning(f"Dose not calculated for ~{round(100*volume_missing/roi['volume'],2)}% of structure {roi['ROI']} ('{struct_names[roi['ROI']]}')")
            else:
                volume_missing = 0
            dose_mean = np.mean(structure_doses)
            integral_dose = col_spacing*row_spacing*frame_spacing*0.001*structure_doses.sum()
        doses, vols = np.unique(structure_doses,return_counts=True)
        # Convert to cumulative dose-volume histogram with relative volume
        vols = np.concatenate(([100],100*np.flip(np.cumsum(np.flip(vols)))/len(structure_doses),[0]))
        doses = np.concatenate(([dose_min],doses,[dose_max]))
        # Calculate dose-volume histogram using 'relative' dose (using fixed 0.5% intervals between min and max doses)
        dvh = _DVH_interp(dvh_dose=doses, dvh_vol=vols, doses=[dose_min+i*0.005*(dose_max-dose_min) for i in range(1,200)])
        logger.info(f"Completed DVH calculation for ROI: {roi['ROI']} ('{struct_names[roi['ROI']]}') [Dmin: {dose_min}, Dmean: {dose_mean}, Dmax: {dose_max})]")
        return_data.append({'ROI': roi['ROI'], 
                            'type': roi['type'], 
                            'dose_min': dose_min, 
                            'dose_mean': dose_mean, 
                            'dose_max': dose_max, 
                            'dose_integral': integral_dose, 
                            'incomplete_dose_volume': volume_missing, 
                            'volume': roi['volume'], 
                            'surface_area': roi['surface_area'], 
                            'dvh': dvh})
    #    #(3004,000A) = Dose Summation Type [Required (1), code string]
    #    dose_sum_type = get_tag(dose,0x3004,0x000A)
    #    #(3004,0005) = Spatial Transform of Dose [Optional (3), code string]
    #    spatial_transform = get_tag(dose,0x3004,0x0005,'NONE')
    # TO DO: incorporate (0070,0404) Referenced Spatial Registration Sequence and (3004,0005) Spatial Transform of Dose!!
    # TO DO: possibly handle multiple dose files if input to sum dose across them (?)
    return(return_data)

def _rois_from_rtstruct(structure, rois=[], res=[0.25,0.25,0.25], force=False, smoothing=0.5, interpolate=2):
    """
    Extracts and processes 3D Regions of Interest (ROIs) from a DICOM RT structure set, considering specified resolution and ROIs.

    This internal function iterates over the ROI Contour Sequence in the provided DICOM RT structure. 
    It processes each ROI based on the provided list of ROIs and resolution parameters. The function supports skipping 
    non-anatomic volumes or technical volumes, depending on the 'force' parameter. For each valid ROI, it extracts the contour 
    data and constructs a 3D mesh representation. It calculates and stores the volume and surface area for each ROI. The function 
    also handles various error cases, logging information about skipped ROIs and exceptions encountered during mesh construction.

    Args:
        structure (pydicom Dataset): The DICOM RT structure dataset, containing structure data.
        rois (list of str, optional): A list of Region of Interest (ROI) names to be processed. Defaults to an empty list, indicating all ROIs.
        res (list of float, optional): A 3-element list specifying the resolution [x, y, z] for mesh construction. Defaults to [0.5, 0.5, 0.5].
        force (bool, optional): If True, processes all ROIs regardless of their type. If False (default), skips non-anatomic or technical ROIs.

    Returns:
        list of dict: A list of dictionaries, each representing an extracted ROI. Each dictionary contains the ROI number, type, mesh object, 
                        volume in cubic centimeters, and surface area in square millimeters.
    """
    return_data = []
    #(0002,0010) = Transfer Syntax UID [Required (1), unique identifier]
    explicit_VR = _get_tag(structure.file_meta,0x0002,0x0010,'') in ['1.2.840.10008.1.2.1','1.2.840.10008.1.2.1.98','1.2.840.10008.1.2.1.99','1.2.840.10008.1.2.2']
    #(3006,0039) = ROI Contour Sequence [Required (1), sequence]
    struct_names = {_get_tag(s,0x3006,0x0022): _get_tag(s,0x3006,0x0026) for s in _get_tag(structure,0x3006,0x0020,[])}
    for roi in _get_tag(structure,0x3006,0x0039,[]):
        #(3006,0084) = Referenced ROI Number [Required (1), integer string]
        ROI = _get_tag(roi,0x3006,0x0084)
        if rois and ROI not in rois:
            logger.info(f"Skipping 3D extraction of ROI: {ROI} ('{struct_names[ROI]}') [missing from expected ROI list]")
            continue
        # Skip processing technical / non-anatomic volumes
        #(3006,0080) = RT ROI Observations Sequence [Required (1), sequence]
        #(3006,0084) = ROI Number [Required (1), integer string]
        #(3006,00A4) = RT ROI Interpreted Type [Required (2), code string] 
        ROI_type = [_get_tag(r,0x3006,0x00A4) for r in _get_tag(structure,0x3006,0x0080,[]) if _get_tag(r,0x3006,0x0084) == ROI][0]
        if not force and ROI_type in ['ISOCENTER','MARKER','SUPPORT','DOSE_MEASUREMENT','FIXATION','CONTROL',
                'BRACHY_ACCESSORY','BRACHY_CHANNEL','BRACHY_CHNL_SHLD','BRACHY_SRC_APP','BOLUS']:
            logger.info(f"Skipping 3D extraction of ROI: {ROI} ('{struct_names[ROI]}', type: {ROI_type})")
            continue
        #(3006,0040) = Contour Sequence [Optional (3), sequence]
        contours = _get_tag(roi,0x3006,0x0040,[])
        if not contours:
            logger.info(f"Skipping 3D extraction of ROI: {ROI} ('{struct_names[ROI]}') [missing contour data]")
            continue
        logger.info(f"Extracting ROI: {ROI} ('{struct_names[ROI]}', type: {ROI_type})...")
        # Extract and post-process structure mesh from DICOM-RT contour data
        try:
            if ROI_type == 'EXTERNAL':
                structure_mesh = _mesh_from_contours(contours,res=max(res,[1,1,1]),explicit_VR=explicit_VR,smoothing=smoothing,interpolate=0)
            else:
                structure_mesh = _mesh_from_contours(contours,res=res,explicit_VR=explicit_VR,smoothing=smoothing,interpolate=interpolate)
        except Exception as e:
            logger.warning(e)
            continue
        mass_properties = vtk.vtkMassProperties()
        mass_properties.SetInputConnection(structure_mesh.GetOutputPort())
        volume_cc = mass_properties.GetVolume()/1000 # mm^3 -> cm^3 (cc)
        surface_area_mm2 = mass_properties.GetSurfaceArea() # mm^2
        return_data.append({'ROI': ROI, 'type': ROI_type, 'mesh': structure_mesh, 'volume': volume_cc, 'surface_area': surface_area_mm2})
    return return_data

def _mesh_from_contours(contours, res=[0.25,0.25,0.25], explicit_VR=False, precision=10, smoothing=0.5, interpolate=2):
    """
    Convert DICOM-RT contours to a VTK mesh (polydata).
    """
    smoothing = np.clip(smoothing,0,1)
    interpolate = np.clip(interpolate,1,10).astype(int)
    padding = 1
    contour_dict = {}
    sweep_lims = {}
    contour_id = 0
    mini_offset = 0.1**(precision+2)
    res = copy.copy(res)
    initial = True
    extents = np.array([float('inf'), float('-inf'), float('inf'), float('-inf'), float('inf'), float('-inf')], dtype=np.float64)
    lims = np.array([float('inf'), float('-inf'), float('inf'), float('-inf'), float('inf'), float('-inf')], dtype=np.float64)
    total_pts = 0
    for contour in contours:
        #(3006,0016) = Contour Image Sequence [Optional (3), sequence]
        #(0008,1155) = Referenced SOP Instance UID [Required (1), unique identifier]
        # TO DO: ???? FIND CORRESPONDING IMAGE AND READING FRAME --> get the image orientation from there along with ?
        #(3006,0042) = Contour Geometric Type [Required (1), code string]
        contour_type = _get_tag(contour,0x3006,0x0042)
        if contour_type != 'CLOSED_PLANAR': # contour_type in ['OPEN_PLANAR','OPEN_NONPLANAR','POINT']:
            continue
        #(3006,0046) = Number of Contour Points [Required (1), integer string]
        n = int(_get_tag(contour,0x3006,0x0046,0))
        total_pts += n
        if not n or n < 3:
            continue
        #(3006,0050) = Contour Data [Required (1), decimal string]
        # Sequence of (x,y,z) triplets defining a contour in the Patient-Based Coordinate System described in DICOM Section C.7.6.2.1.1 (mm).
        pts = _get_tag(contour,0x3006,0x0050,[])
        # (Potentially) improper contour data encoding if Explicit VR transfer syntax and VL > 65534 bytes
        if explicit_VR and sum(p.__sizeof__() for p in pts) > 65534:
            raise Exception(f'Skipping 3D extraction of ROI: {ROI} (improper contour data encoding)')
        contour_pts = np.array(pts, dtype=np.float64).reshape(-1, 3)
        if initial:
            a, b, c, d = _compute_plane_from_points(contour_pts)
            if d is None:
                continue
            # Set sweep=0, 1, or 2 if the plane is most orthogonal to the z, x, or y axis, respectively
            sweep = (np.argmax(np.array([abs(a), abs(b), abs(c)]))+1) % 3
            indices = {0: [0, 1], 1: [1, 2], 2: [2, 0]}[sweep]
            alt = (sweep+1) % 3
            other = 3-sweep-alt
            res_sweep = res[sweep]
            res_alt = res[alt]
            res_other = res[other]
            res_x, res_y, res_z = res
            if [a,b,c][other] != 1:
                raise Exception(f'Skipping DVH calculation for ROI: {ROI} (structure not defined on an x, y, or z plane)')
            normal = np.array([a, b, c])
            normal /= np.linalg.norm(normal)
            initial = False
        else:
            # Ensure parallel planes
            a, b, c, d = _compute_plane_from_points(contour_pts, a, b, c)
            if d is None:
                continue
            if not (0.999 <= abs(np.dot(normal, np.array([a, b, c])/np.linalg.norm(np.array([a, b, c])))) <= 1.001):
                raise Exception(f'Skipping DVH calculation for ROI: {ROI} (structure defined on non-parallel planes)')
        if np.array_equal(contour_pts[0], contour_pts[-1]):
            contour_pts = contour_pts[:-1]
        plane_key = (round(a,2), round(b,2), round(c,2), round(d,2))
        if plane_key not in contour_dict:
            contour_dict[plane_key] = {'plane': (a, b, c, d), 'contour_pts': [contour_pts]}
        else:
            contour_dict[plane_key]['contour_pts'].append(contour_pts)
        min_vals = np.amin(contour_pts, axis=0)
        max_vals = np.amax(contour_pts, axis=0)
        lims = np.array([np.minimum(lims[0],min_vals[0]),np.maximum(lims[1],max_vals[0]),
                np.minimum(lims[2],min_vals[1]),np.maximum(lims[3],max_vals[1]),
                np.minimum(lims[4],min_vals[2]),np.maximum(lims[5],max_vals[2])], dtype=np.float64)
    extents = np.array([np.round(np.floor(lims[0]/res_x)*res_x,precision), np.round(np.ceil(lims[1]/res_x)*res_x,precision), 
                np.round(np.floor(lims[2]/res_y)*res_y,precision), np.round(np.ceil(lims[3]/res_y)*res_y,precision), 
                np.round(np.floor(lims[4]/res_z)*res_z,precision), np.round(np.ceil(lims[5]/res_z)*res_z,precision)], dtype=np.float64)
    # Sort plane keys into an ordered stack (assumption: contour planes are all parallel)
    plane_keys = sorted(contour_dict.keys(), key=lambda x: x[3])
    if len(plane_keys) == 0:
        raise Exception(f'Skipping DVH calculation for ROI: {ROI} (empty structure)')
    # Fill in any gaps with missing plane keys (assumption: planes are evenly spaced)
    diffs = [plane_keys[i+1][3] - plane_keys[i][3] for i in range(len(plane_keys)-1)]
    d_diff = set(diffs)
    if len(d_diff) > 1:
        d_diff = min(d_diff)
        list(enumerate(diffs))[::-1]
        for d, diff in list(enumerate(diffs))[::-1]:
            if diff > d_diff:
                num_missing = int(diff / d_diff) - 1
                for j in range(num_missing):
                    plane_keys.insert(d+j+1, (plane_keys[d][0],plane_keys[d][1],plane_keys[d][2],plane_keys[d][3]+d_diff*(j+1)))
    n_planes = len(plane_keys)
    sweep_ext = extents[(sweep*2):(sweep*2+2)]
    a_ext = extents[(alt*2):(alt*2+2)]
    o_ext = extents[(other*2):(other*2+2)]
    n_sweeps = int(round(sweep_ext[1]-sweep_ext[0],precision)/res_sweep)
    n_a = int(round(a_ext[1]-a_ext[0],10)/res_alt)
    # Ensure minimum resolution to improve accuracy in processing small structures
    while n_sweeps < 100:
        res_sweep *= 0.5
        res[sweep] *= 0.5
        n_sweeps = int(round(sweep_ext[1]-sweep_ext[0],precision)/res_sweep)
    while n_a < 100:
        res_alt *= 0.5
        res[alt] *= 0.5      
        n_a = int(round(a_ext[1]-a_ext[0],10)/res_alt)
    # Establish list of midpoint coordinates for all x, y cells in voxel_data
    # Establish list of midpoint coordinates for all x, y cells in voxel_data
    sweep_vals = np.round(np.linspace(sweep_ext[0]+0.5*res_sweep, sweep_ext[1]-0.5*res_sweep, n_sweeps), precision)
    alt_vals = np.round(np.linspace(a_ext[0]+0.5*res_alt, a_ext[1]-0.5*res_alt, n_a), precision)
    alt_grid, sweep_grid = np.meshgrid(alt_vals, sweep_vals)
    sweep_flat = sweep_grid.ravel()
    alt_flat = alt_grid.ravel()
    pln_points = np.column_stack((sweep_flat, alt_flat))
    res[other] = (o_ext[1]-o_ext[0])/n_planes
    np_voxels = np.full([n_sweeps+2*padding, # 'sweep' dimension (if sweep=0: x, sweep=1: y, sweep=2: z)
                    n_a+2*padding, # 'alt' dimension (if sweep=0: y, sweep=1: z, sweep=2: x)
                    (n_planes-1)*interpolate+1+2*padding], # 'other' dimension (if sweep=0: z, sweep=1: x, sweep=2: y)
                    res[other]) # fill array with this value
    # Iterate through the ordered plane_keys
    for i, pln in enumerate(plane_keys):
        pln_index = (n_planes-i-1)*interpolate+padding
        if pln not in contour_dict:
            continue
        vertices = []
        codes = []
        for contour in contour_dict[pln]['contour_pts']:
            codes.extend([Path.MOVETO] + [Path.LINETO] * contour.shape[0])
            vertices.append(np.vstack([contour, contour[0]]))
        poly_set = Path(np.concatenate(vertices, axis=0)[:, indices], codes)        
        inside_pts = poly_set.contains_points(pln_points).reshape((n_sweeps, n_a))
        sweep_indices, alt_indices = np.floor((poly_set.vertices-np.array([sweep_ext[0], a_ext[0]]))/np.array([res_sweep, res_alt])).astype(int).T
        valid_indices = (sweep_indices >= 0) & (sweep_indices < inside_pts.shape[0]) & (alt_indices >= 0) & (alt_indices < inside_pts.shape[1])
        sweep_indices = sweep_indices[valid_indices]
        alt_indices = alt_indices[valid_indices]
        inside_pts[sweep_indices, alt_indices] = 1
        np_voxels[padding:padding+n_sweeps,padding:padding+n_a,pln_index] = np.where(inside_pts, -ndimage.distance_transform_edt(inside_pts), ndimage.distance_transform_edt(~inside_pts))
        #plt.imshow(np_voxels[:,:,pln_index], cmap='jet', interpolation='nearest')
        #plt.colorbar()
        #plt.contour(np.abs(np_voxels[:,:,pln_index]) < 1.05, colors='black', levels=[0.95,1.05])
        #plt.show()
    voxel_dims = np_voxels.shape
    slice_offset = np.full((voxel_dims[0],voxel_dims[1]), -res[other])
    if interpolate > 1:
        res[other] /= interpolate
        for pln_index in range(padding,np_voxels.shape[2]-padding-interpolate+1,interpolate):
            SDF_a = np_voxels[:,:,pln_index]
            SDF_b = np_voxels[:,:,pln_index+interpolate]
            for j in range(1,interpolate):
                t = j/interpolate
                np_voxels[:,:,pln_index+j] = SDF_b*t+SDF_a*(1-t)
    np_voxels = np.transpose(np_voxels,(other,alt,sweep))
    voxel_data = vtk.vtkImageData()
    voxel_data.SetSpacing(res)
    # TO DO: NEED TO VERIFY ORIGINS RECONCILE W/ DOSE GRID COORDS,etc.
    voxel_data.SetOrigin(extents[0]-padding*res[0], extents[2]-padding*res[1], extents[4]-padding*res[2])
    voxel_data.SetDimensions(np_voxels.shape[::-1])
    voxel_data.GetPointData().SetScalars(vtk_np.numpy_to_vtk(num_array=np_voxels.ravel(), deep=True, array_type=vtk.VTK_FLOAT))    
    # Use Flying Edges (3D) algorithm to identify structure surface
    surface = vtk.vtkFlyingEdges3D()
    surface.SetInputData(voxel_data)
    surface.SetValue(0, 0.0)
    surface.ComputeNormalsOff()
    surface.ComputeGradientsOff()
    surface.ComputeScalarsOff()
    surface.Update()
    # Reduce surface mesh size for larger surfaces (>25,000 triangles) via decimation to improve downstream smoothing and stencil computation efficiency
    n_triangles = surface.GetOutput().GetNumberOfCells()
    if n_triangles > 25000:
        decimator = vtk.vtkDecimatePro()
        decimator.SetInputConnection(surface.GetOutputPort())
        # Increase decimation ratio logarithmically from 0 to 0.95 as surface mesh size increases to >5 million triangles
        decimator.SetTargetReduction(0.95*np.log10(min(n_triangles,5000000)/25000)/np.log10(200))
        decimator.PreserveTopologyOn()
        decimator.SetFeatureAngle(60)
        decimator.SetMaximumError(1)
        decimator.SplittingOff()
        decimator.Update()
    else:
        decimator = surface
    # Smooth surface mesh using Windowed Sinc PolyData filter
    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputConnection(decimator.GetOutputPort())
    # Transform smoothing factor to NumberOfIterations for WindowedSinc, according to following scale:
    #   smoothing=0.0  ->  surface_nets_iters=0   (almost no smoothing)
    #   smoothing=0.2  ->  surface_nets_iters=2   (little smoothing)
    #   smoothing=0.5  ->  surface_nets_iters=8   (average smoothing)
    #   smoothing=0.7  ->  surface_nets_iters=14  (strong smoothing)
    #   smoothing=1.0  ->  surface_nets_iters=24  (very strong smoothing)
    smoother.SetNumberOfIterations(round(15*smoothing**2+9*smoothing))
    smoother.BoundarySmoothingOff()
    smoother.FeatureEdgeSmoothingOff()
    # Transform smoothing factor to PassBand for WindowedSinc, according to following scale:
    #   smoothing=0.0  ->  pass_band=1.0    (almost no smoothing)
    #   smoothing=0.25 ->  pass_band=0.1    (average smoothing)
    #   smoothing=0.5  ->  pass_band=0.01   (more smoothing)
    #   smoothing=0.75 ->  pass_band=0.001  (strong smoothing)
    #   smoothing=1.0  ->  pass_band=0.0001 (very strong smoothing)
    smoother.SetPassBand(pow(10.0,-4*smoothing))
    smoother.NonManifoldSmoothingOn()
    smoother.NormalizeCoordinatesOff()
    smoother.Update()
    # Return original surface if Windowed Sinc smoothing fails
    if any(abs(b) > 1e20 for b in smoother.GetOutput().GetBounds()):
        return(decimator)
    return(smoother)

def _stencil_from_mesh(meshdata_port, dose_image, dose_orientation=None):
    """
    Generate a stencil from a VTK mesh.
    """
    extent = dose_image.GetExtent()
    spacing = dose_image.GetSpacing()
    origin = dose_image.GetOrigin()
    # Transform mesh coordinates to dose grid coordinates using inverse transform
    if dose_orientation:
        transform = vtk.vtkTransform()
        transform.SetMatrix(dose_orientation)
        transform.Inverse()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetInputConnection(meshdata_port)
        transform_filter.SetTransform(transform)
        transform_filter.Update()
        meshdata_port = transform_filter.GetOutputPort()
    stencil = vtk.vtkPolyDataToImageStencil()
    stencil.SetInputConnection(meshdata_port)
    stencil.SetOutputOrigin(origin)
    stencil.SetOutputSpacing(spacing)
    stencil.SetOutputWholeExtent(extent)
    stencil.Update()
    return stencil

def _compute_plane_from_points(points, a=None, b=None, c=None):
    """
    Compute the plane equation Ax + By + Cz + D = 0 given at least 3 non-collinear points in the space.
    If the points are collinear or not enough points are provided, return None.
    Parameters:
        points (list): A list of 3D points [(x1, y1, z1), (x2, y2, z2), ...].

    Returns:
        tuple: Coefficients (a, b, c, d) for the plane equation, or (None, None, None, None) if the points are collinear.
    """
    if len(points) < 1:
        return a, b, c, None
    if a is not None and b is not None and c is not None:
        d = -a*points[0][0]-b*points[0][1]-c*points[0][2]
    else:
        d = None
    if len(points) < 3:
        return a, b, c, d
    points = np.array(points)
    for i in range(min(len(points) - 2,10)):
        p1 = points[i]
        p2 = points[i + 1]
        p3 = points[i + 2]
        # Compute normal of the plane
        normal = np.cross(np.subtract(p2, p1), np.subtract(p3, p1))
        # Check if the points are collinear
        if np.linalg.norm(normal) != 0:
            normal = normal / np.linalg.norm(normal)  # normalize
            # Adjust the normal to ensure z-component is positive
            if normal[2] < 0:
                normal = -normal
            # Compute the D coefficient
            D = -np.dot(normal, p1)
            return tuple(normal) + (D,)
    if d:
        return a, b, c, d
    # If no non-colinear sets of points, compute the best fit plane using PCA
    centroid = points.mean(axis=0)
    centered_data = points - centroid
    cov_matrix = np.cov(centered_data, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    normal = eigenvectors[:, np.argmin(eigenvalues)]
    D = -np.dot(normal, centroid)
    return tuple(normal) + (D,)

def _DVH_interp(dvh_dose, dvh_vol, doses=None, volumes=None, vol_min=0, vol_max=100):
    """Linearly interpolate dose-volume histogram (DVH) data at arbitrary points
    dvh_dose: list of numerical values corresponding to DVH dose points
    dvh_vol: list of numerical values corresponding to DVH volume points
    doses: list of numerical values at which to return interpolated volume data
    volumes: list of numerical values at which to return interpolated dose data
    vol_min: extreme 'left' value to return (avoid extrapolation)
    vol_max: extreme 'right' value to return (avoid extrapolation)
    Return value: List containing one or more interpolated dose values
    """
    if doses:
        return [
            vol_max if not i 
            else vol_min if i == len(dvh_dose) 
            else dvh_vol[i-1]+(d-dvh_dose[i-1])*(dvh_vol[i]-dvh_vol[i-1])/(dvh_dose[i]-dvh_dose[i-1])
            for d in doses
            for i in [bisect_left(dvh_dose, d)]
        ]
    elif volumes:
        dvh_vol = dvh_vol[::-1]
        dvh_dose = dvh_dose[::-1]
        return [
            dvh_dose[0] if not i
            else dvh_dose[-1] if i==len(dvh_vol)
            else dvh_dose[i-1]+(v-dvh_vol[i-1])*(dvh_dose[i]-dvh_dose[i-1])/(dvh_vol[i]-dvh_vol[i-1])
            for v in volumes
            for i in [bisect_left(dvh_vol, v)]
        ]

__all__ = ['dicom_validation', 'dicom_extract', 'match_pt_identity_from_dicom', 'calculate_dvh']