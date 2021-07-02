import numpy as np

ICMW_UWLCM_names = {
  "time" : "position",
  "T" : "T",
  "Qv" : "Qv",
  "RH" : "RH",
  "S_drop" : "S_drop",
  "LWC" : "LWC_gm-3",
  "N_drop" : "N_drop",
  "N_aerosol" : "N_aerosol",
  "N_removal" : "N_removal",
  "disp_r" : "disp_r",
  "r_mean1" : "r_mean1",
  "Sigma2_S" : "Sigma2_S",
  "Sigma2_T" : "Sigma2_T",
  "Sigma2_Qv" : "Sigma2_Qv",
  "Sigma2_S_drop" : "Sigma2_S_drop",
  "epsilon" : "epsilon",
  "TKE" : "tot_tke",
  "H_flux_t" : "H_flux_t",
  "H_flux_b" : "H_flux_b",
  "qv_flux_t" : "qv_flux_t",
  "qv_flux_b" : "qv_flux_b"
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

#example use:
for key in ICMW_UWLCM_names:
  print key
  print read_ICMW_var(open("/home/piotr/praca/ICMW_Pi_chamber_LES/wyniki/series/output_pichamber_3D_sgs_turbadve0_source1e5SD1@domainSstp70_SideRH80_Coal0_out_lgrngn_pi_chamber_icmw_series.dat","r"), key)
