import os
import sys
#from skimage import transform
from pydicom import dcmread
import numpy as np
import argparse
import torch
#from data import create_dataset
import pdb
import re
import pandas as pd
import scipy.io as sio
from scipy.interpolate import interp1d

sys.path.append('.')


def get_torch_tensor(npy_tensor, device):
    out = torch.from_numpy(npy_tensor)
    out.to(device)

    return out

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
    
    file_caseinfo = "Z:\LulinY\Lung-dosimetrics\\2024\python_code\AI_RTP\VCU_Lung_plan_info_10_17_2024.csv"
    case_info = pd.read_csv(file_caseinfo)
    
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
    column_names = ['CaseID', 'PlanID', 'Dx'] + oars2
    df_cases = pd.DataFrame(columns=column_names)
    filename_caseinfo = os.path.join(out_dir, 'Cases_processing_summary.csv')
 
    #fh = open(filename,"w",newline="")

    #for idx, case in enumerate(cases):
    for idx in range(0, len(cases)):
        case = cases[idx]
        if not case.endswith('.npz'):
            continue
 
        case_path = os.path.join(in_dir, case)
        caseID = re.search(r'^VCU_Lung_\d+_\d', case).group()
        patID = caseID[0:-2]
        planID = caseID[-1]
        w311 = case_info['CaseID'].str.contains(patID) 
        w312 = case_info['PlanID'].astype('string').str.contains(planID)
        w31 = w311 & w312
        if any(w31):
            Dx = case_info.PrescriptionDose[w31].iloc[0]
            Compat = case_info.Compatible_with_DoseRTX[w31].iloc[0]
        else:
            continue
        
        if np.isnan(Compat) or (not bool(Compat)):
            continue
 
        fh = open(filename_caseinfo,"w",newline="")        
        row_data1 = {'CaseID':caseID, 'PlanInd':planID, 'Dx': Dx}
        
        try:
            print('Processing case {}: {} of {} ...'.format(case, idx+1, len(cases)))
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
            
            num_ptv = np.sum(ptv)
            w20 = dict_np['DOSE']
            dose_copy = w20.copy()
            dose_copy *= ptv
            #return ptv, oar
            #print(ptv)
            #print(dose_copy)
            sum = np.sum(dose_copy)
            scale_factor = (60 * num_ptv) / sum    ##debug
            w20 *= scale_factor    ##debug

            
            n_bins = 376                               ## for his_sig with 376 bins
            bins = np.linspace(0, 75, n_bins)
            bins_scaled = bins*scale_factor
            w13 = dict_np['HIST_SIG']
            w14 = np.delete(w13,[3,5],axis=1)
            f1 = interp1d(bins_scaled, w14,kind='linear',axis=0, fill_value=0, bounds_error=False)
            w141 = f1(bins)
            
            bins = dict_np['BINS']
            bins_scaled = bins*scale_factor
            w15 = dict_np['HIST_CLINICAL']
            w16 = np.delete(w15,[3,5],axis = 1)
            f2 = interp1d(bins_scaled, w16,kind='linear',axis=0, fill_value=0, bounds_error=False)
            w161 = f2(bins)

# =============================================================================
#             if Dx>0.001:
#                 w21 = w20/Dx
#             else:
#                 continue
#             
# =============================================================================

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
            df_cases.loc[len(df_cases)] = row_data  
            
            filename_npz = os.path.join(out_dir, caseID+'_reduced.npz')

            np.savez(filename_npz, CT=dict_np['CT'], DOSE=w20, OAR=w12, PTV=ptv, HIST_SIG=w141, HIST_CLINICAL = w161, BINS=dict_np['BINS'])
            filename_mat = os.path.join(out_dir, caseID+'_reduced.mat')
            sio.savemat(filename_mat, {'DOSE': w20, 'OAR': w12, 'HIST_REAL': w161, 'BINS':dict_np['BINS']})     
            # Export the DataFrame as a CSV spreadsheet
            df_cases.to_csv(fh)
            fh.flush()
            fh.close()
        except:
            print('Processing case {} failed...'.format(caseID))
            fh.close()
    #fh.close()    
    # debug
    return df_cases, w16
## DEbug
#out_dir = "Z:\LulinY\Lung-dosimetrics\2024\python_code\DoseRTX\VCU_Lung_2024_processed_data"
in_dir = "C:\Lulin-home\KBP-lung\CE project\AI_RTP\VCU_Lung_2024_processed_data"
in_dir = "G:/My Drive/AI_RTP/VCU_Lung_preprocess"
in_dir = "R:\LulinY\Processed_NPZ"
out_dir = "C:\Lulin-home\KBP-lung\CE project\AI_RTP\VCU_Lung_2024_dataset\test"
out_dir = "G:/My Drive/AI_RTP/VCU_Lung_preprocess/test"
out_dir = "R:\LulinY\Processed_NPZ_reduced"
w2, hist = batch_process_cases(in_dir,out_dir)
#fh.close()