import numpy as np

ICMW_UWLCM_names = {
  "time" : "position",
  "T" : "T_nowall",
  "Qv" : "Qv_nowall",
  "RH" : "RH_nowall",
  "LWC" : "LWC_nowall",
  "LWC_gm-3" : "LWC_gm-3_nowall",
  "N_drop" : "N_drop_nowall",
  "N_aerosol" : "N_aerosol_nowall",
  "N_removal" : "N_removal",
  "disp_r" : "disp_r_nowall",
  "r_mean1" : "r_mean1_nowall",
  "r_mean2" : "r_mean2_nowall",
  "Sigma2_S" : "Sigma2_S_nowall",
  "Sigma2_T" : "Sigma2_T_nowall",
  "Sigma2_Qv" : "Sigma2_Qv_nowall",
  "epsilon" : "epsilon_nowall"
}

def read_UWLCM_array(file_obj):
  arr_name = file_obj.readline()
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr), arr_name

def read_UWLCM_var(file_obj, var_name):
  file_obj.seek(0)
  while True:
    arr, name = read_UWLCM_array(file_obj)
    if(str(name).strip() == str(var_name).strip()):
      break
  return arr

def read_ICMW_var(file_obj, ICMW_var_name):
  return read_UWLCM_var(file_obj, ICMW_UWLCM_names[ICMW_var_name])

for key in ICMW_UWLCM_names:
  print key
  print read_ICMW_var(open("/home/piotr/praca/ICMW_Pi_chamber_LES/wyniki/series_and_profs/output_pichamber_3D_sgs_turbadve0_source5e9SD1000@center_Coal0_out_lgrngn_pi_chamber_icmw_series.dat","r"), key)
