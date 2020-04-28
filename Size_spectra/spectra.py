import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import OrderedDict

plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
plot_dry = int(argv[4])
outfile = argv[5]
#from_lvl = int(argv[4])
#to_lvl = int(argv[5])

if plot_dry == True:
  size_data = {"rw" : 29, "rd" : 21}
else:
  size_data = {"rw" : 29}

#hardcoded bin edges, need to match UWLCM
left_edges = {"rw": np.zeros(30), "rd": np.zeros(22)}
bin_centers = {"rw": np.zeros(29), "rd": np.zeros(21)}
bin_width = {"rw": np.zeros(29), "rd": np.zeros(21)}

for i in np.arange(30):
  left_edges["rw"][i] = 10**(-3 + i * .2) * 1e-6 # [m]

for i in np.arange(29):
  bin_centers["rw"][i] = 0.5 * (left_edges["rw"][i] + left_edges["rw"][i+1])
  bin_width["rw"][i] = (left_edges["rw"][i+1] - left_edges["rw"][i])

for i in np.arange(22):
  left_edges["rd"][i] = 10**(-3 + i * .2) * 1e-6 # [m]
for i in np.arange(21):
  bin_centers["rd"][i] = 0.5 * (left_edges["rd"][i] + left_edges["rd"][i+1])
  bin_width["rd"][i] = (left_edges["rd"][i+1] - left_edges["rd"][i])


print left_edges
print bin_centers
print bin_width

data_names = {}
for rwrd in size_data:
  bin_no = np.arange(0,size_data[rwrd])
  data_names[rwrd] = []
  for no in bin_no:
    data_names[rwrd] = np.append(data_names[rwrd], rwrd + "_rng" + str(no).zfill(3) + "_mom0")

print data_names

layer_thickness = 10
cloud_thresh = 1e-8

directories = argv[6:len(argv):2]
labels = argv[7:len(argv):2]
print directories, labels

levels = ["ground", "cloud_base"]
#levels = ["all", "pi_chamber_measurement_location"]

if plot_dry == True:
  all_data_names = np.append(data_names["rw"], data_names["rd"])
else:
  all_data_names = data_names["rw"]

for lvl in levels:
  # range of radii plotted
  # full range
  #if plot_dry == True:
  #  r_min = min(left_edges["rw"][0], left_edges["rd"][0])
  #  r_max = max(left_edges["rw"][29], left_edges["rd"][21])
  #else:
  #  r_min = left_edges["rw"][0]
  #  r_max = left_edges["rw"][29]
  # range adapted for nonzero values, init with some crazy values that will make any values found smaller/larger
  r_min = 1000
  r_max = -1000

  total_arr = OrderedDict()
  for data in all_data_names:
    total_arr[data] = OrderedDict()
  
  plot_labels = OrderedDict()
  tot_cloud_base_lvl = OrderedDict()
  for lab in labels:
    tot_cloud_base_lvl[lab] = np.zeros(0)
  
  # read in nx, ny, nz
  for directory, lab in zip(directories, labels):
    rhod = h5py.File(directory + "/const.h5", "r")["G"][:,:,:]
    w3d = h5py.File(directory + "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:,:]
    nx, ny, nz = w3d.shape
    plot_labels[lab] = lab
  #  Exy_avg = OrderedDict()
  #  for data in data_names:
  #    Exy_avg[data] = np.zeros(((nx+1)/2))
  
    for data in all_data_names:
      total_arr[data][lab] = np.zeros(0)
  
      for t in range(time_start, time_end+1, outfreq):
        filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
        #print filename

        if lvl == "cloud_base":
          # find cloud base
          # based on cloud rw
          w3d = h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:] # * 4. / 3. * 3.1416 * 1e3 
          cloud_base_lvl = np.argmax(np.average(w3d, axis=(0,1)) > cloud_thresh)
          # based on RH
    #      w3d = h5py.File(filename, "r")["RH"][:,:,:] # * 4. / 3. * 3.1416 * 1e3 
    #      cloud_base_lvl = np.argmax(np.average(w3d, axis=(0,1)) > .99)
          tot_cloud_base_lvl[lab] = np.append(tot_cloud_base_lvl[lab], cloud_base_lvl) # done for each data, but we dont care - wont affect average
          print 'cloud base lvl = ', cloud_base_lvl
          total_arr[data][lab] = np.append(total_arr[data][lab], (h5py.File(filename, "r")[data]*rhod)[:,:,cloud_base_lvl-layer_thickness : cloud_base_lvl])
        if lvl == "ground":
          total_arr[data][lab] = np.append(total_arr[data][lab], (h5py.File(filename, "r")[data]*rhod)[:,:, 0 : layer_thickness ])
        if lvl == "all":
          total_arr[data][lab] = np.append(total_arr[data][lab], (h5py.File(filename, "r")[data]*rhod)[:,:,:])
        if lvl == "pi_chamber_measurement_location":
          total_arr[data][lab] = np.append(total_arr[data][lab], (h5py.File(filename, "r")[data]*rhod)[32:34,5:6,4:5]) # roughly + assuming dx=0.03125 m
  
  #    hists[lab] = np.hist(total_arr, bins=100)
  #    _ = plt.hist(total_arr, bins='auto')
  
  for lab in labels:
  #  print  np.average(total_arr[lab])
    plot_labels[lab] = plot_labels[lab] + ' <cloud base lvl> = {:.2f}'.format(np.average(tot_cloud_base_lvl[lab]))
    
  plt.rcParams.update({'font.size': 30})
  plt.figure(figsize=(40,40))
  #_ = plt.hist(total_arr["rain_rw_mom3"].values(), bins='auto', label=plot_labels.values(), density=True)

  avg_conc = {}
  avg_conc_arr = {}
  for rwrd in size_data:
    avg_conc_arr[rwrd] = {}
    for lab in labels:
      avg_conc_arr[rwrd][lab] = np.zeros(size_data[rwrd])
    for name,it in zip(data_names[rwrd], np.arange(0, size_data[rwrd])):
      avg_conc[name] = {}
      for lab in labels:
        avg_conc[name][lab] = np.average(total_arr[name][lab])
        avg_conc_arr[rwrd][lab][it] =  avg_conc[name][lab]

    for lab in labels:
      #print avg_conc_arr[rwrd][lab]
      nonzero_indices = [index for index, item in enumerate(avg_conc_arr[rwrd][lab]) if item != 0]
      print nonzero_indices
      first_nonzero_idx = nonzero_indices[0]
      last_nonzero_idx = nonzero_indices[-1]
      r_min = min(r_min, left_edges[rwrd][first_nonzero_idx])
      r_max = max(r_max, left_edges[rwrd][last_nonzero_idx+1])
      print 'r_min = ', r_min, ' r_max = ', r_max
      plt.step(bin_centers[rwrd] * 1e6 * 2, avg_conc_arr[rwrd][lab] / bin_width[rwrd] / 1e12 / 2, where='mid', label=rwrd + '_' + lab, linewidth=6) # *1e6 to have microns on x, / 1e12 to adjust for width in microns and to have concentration per cm^3; *2 and /2 to get diameters

#    data = total_arr["rain_rw_mom3"].values()

  plt.xlabel('diameter [um]')
  plt.ylabel('PDF of concentration [cm^{-3} / um]')
  plt.xlim(r_min*1e6*2, r_max*1e6*2)
 # plt.xscale('log')
  plt.legend()
  plt.yscale('log')
  plt.savefig(outfile + '_size_spectra_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
  
