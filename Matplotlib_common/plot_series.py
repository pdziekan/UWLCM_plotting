import numpy as np
from sys import argv
from math import floor

from latex_labels import var_labels
from read_UWLCM_arrays import *

def plot_series(var_list, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict, ylimdict, show_bin=False, suffix='', xlabel='', ylabeldict=var_labels, file_names=[], file_labels=[], linewidth=1):
  # if file names are not defined, read them and labels from command line
  if len(file_names)==0:
    file_no = np.arange(1, len(argv)-1 , 2)
    for no in file_no:
      file_names.append(argv[no] + suffix)
      file_labels.append(argv[no+1])

  for var in var_list:
    label_counter=0
    for file_name in file_names:
      print(file_name, var)
      series_file = open(file_name, "r")
      my_times = read_my_var(series_file, "position")
      my_res = read_my_var(series_file, var)
      
      # rescale time to hours
      my_times = my_times / 3600.

      # rescale autoconv. and accr. rates to g/(m^3 * day)
      if var == "cl_accr25_rico" or var == "cl_acnv25_rico"  or var == "cl_acnv25_dycoms" or var == "cl_accr25_dycoms":
        my_res *= 3600 * 24
      
      series_file.close()
  
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
      colorList = ['red', 'blue', 'green']

      # x label only on he lowest row
      xlabel_used = xlabel
      if plot_iter < nploty:
        xlabel_used = ''

      plot_my_array(axarr, plot_iter, my_times, my_res, nploty, xlabel=xlabel_used, ylabel=ylabeldict[var], varlabel=file_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xscale=xscaledict[var], yscale=yscaledict[var], xlim=xlimdict[var], ylim=ylimdict[var], linewidth=linewidth)#, color = colorList[int(floor(label_counter / len(dashList)))])
      label_counter+=1
    plot_iter = plot_iter + 1
  return plot_iter
