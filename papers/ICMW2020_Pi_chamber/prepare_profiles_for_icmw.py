import numpy as np
import os

ICMW_UWLCM_names = {
  "time" : "position",
  "T" : "T",
  "Qv" : "rv",
  "RH" : "RH"
}

inj_rate_to_N = {
  "0.25e9" : "13",
  "0.5e9"  : "35",
  "0.75e9" : "60",
  "1.2e9"  : "120",
  "1.5e9"  : "180",
  "2.5e9"  : "400",
  "5e9"    : "1200",
  "6e9"    : "1800",
  "7.5e9"    : "2500"
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

time_points = np.arange(0,3601,30)

for inj_rate in inj_rate_to_N:
  for var in ["T", "RH", "Qv"]:
    out_file_name = "/scratch/plgdziekan/LES_N"+inj_rate_to_N[inj_rate]+"_2D_"+var
    if os.path.exists(out_file_name):
      os.remove(out_file_name)
    out_file=open(out_file_name,'a')
    for time in time_points:
      arr = read_UWLCM_var(open("/scratch/plgdziekan/output_pichamber_3D_sgs_turbadve0_source"+inj_rate+"SD1000@center_Coal0_out_lgrngn_pi_chamber_icmw_profiles_"+str(time)+"_"+str(time)+".dat","r"), ICMW_UWLCM_names[var])
      np.savetxt(out_file, arr[5:28], newline=" ", delimiter=" ", fmt="%1.5f")
      out_file.write("\n")
