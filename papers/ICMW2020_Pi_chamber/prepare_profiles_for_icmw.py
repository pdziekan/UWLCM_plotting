import numpy as np
import os
from sys import argv

ICMW_UWLCM_names = {
  "time" : "position",
  "T" : "T",
  "Qv" : "rv",
  "RH" : "RH"
}

inj_rate_to_N = {
  "1e5" : "4",
  "3e5" : "18",
  "1e6" : "120",
  "3e6" : "700",
}

inj_rate_SstpSrc = {
  "1e5" : "70",
  "3e5" : "70",
  "1e6" : "100",
  "3e6" : "100",
}

def read_UWLCM_array(file_obj):
  arr_name = file_obj.readline()
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = [float(x) for x in line]
  return np.array(arr), arr_name

def read_UWLCM_var(file_obj, var_name):
  file_obj.seek(0)
  while True:
    arr, name = read_UWLCM_array(file_obj)
    if(str(name).strip() == str(var_name).strip()):
      break
  return arr

time_points = np.arange(0,5401,30)

assert(len(argv)==2)

data_dir = argv[1]

for inj_rate in inj_rate_to_N:
  for var in ["T", "RH", "Qv"]:
    out_file_name = data_dir+"/LES_N"+inj_rate_to_N[inj_rate]+"_2D_"+var
    if os.path.exists(out_file_name):
      os.remove(out_file_name)
    out_file=open(out_file_name,'a')
    for time in time_points:
      file_name = data_dir+"/output_pichamber_3D_sgs_turbadve0_source"+inj_rate+"SD1@domainSstp"+inj_rate_SstpSrc[inj_rate]+"_SideRH80_Coal0_out_lgrngn_pi_chamber_icmw_profiles_"+str(time)+"_"+str(time)+".dat"
      arr = read_UWLCM_var(open(file_name,"r"), ICMW_UWLCM_names[var])
      np.savetxt(out_file, arr[:], newline=" ", delimiter=" ", fmt="%1.5f")
      out_file.write("\n")
