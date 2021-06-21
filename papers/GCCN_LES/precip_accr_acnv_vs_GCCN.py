from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
from read_UWLCM_arrays import read_my_var
from latex_labels import labeldict

# activate latex text rendering
rc('text', usetex=True)

series_from_it = int(sys.argv[1])
series_to_it = int(sys.argv[2])
profs_from_it = int(sys.argv[3])
profs_to_it = int(sys.argv[4])
qlimit = float(sys.argv[5])

#varlabels = ["{\it ScNc30}", "{\it ScNc40\_salt\_CCN}", "{\it ScNc45}", "{\it ScNc105}"]
#varlabels = ["{\it Cu38}", "{\it Cu60}", "{\it Cu85}"]
varlabels = ["{\it Sc38}", "{\it Sc60}", "{\it Sc115}"]
averaging_period = float(profs_to_it - profs_from_it) / 3600. # period over which series are averaged [h]; NOTE: we assume that series_from(to)_it = profs_from(to)_it / outfreq!

# assumed initial GCCN concentrations
GCCN_conc = [0,0.216,5*0.216,10*0.216]

# init the plot
nplotx = 2 #int(nplots/6 + 0.5)
nploty = 2
fig, axarr = plt.subplots(nploty, nplotx)#), constrained_layout=True )

#prepare a list of output files, assuming the following order: prsitine, standard, polluted, for each 4 results: no GCCN, GCCN ,GCCNx5, GCCNx10
series_file_names = []
profs_file_names = []
file_no = np.arange(6, len(sys.argv)-1 , 1)
print file_no
for no in file_no:
  series_file_names.append(sys.argv[no] + "series.dat")
  profs_file_names.append(sys.argv[no] + "profiles_"+str(profs_from_it)+"_"+str(profs_to_it)+".dat")

#assert(len(series_file_names) == 16)

label_counter=0

for it in np.arange(12):
  print(series_file_names[it])
  print(profs_file_names[it])

  if(it % 4 == 0):
    mean_surf_precip = []
    tot_acc_surf_precip = []
    tot_acc_surf_precip_std_dev = []
    prflux = []
    prflux_std_dev = []
    prfluxDivByClFrac= []
    prfluxDivByClFrac_std_dev = []
    prfluxFromPrfluxVsClhght = []
    prfluxFromPrfluxVsClhght_std_dev = []
    tot_acc_acnv = []
    tot_acc_acnv_std_dev = []
    tot_acc_accr = []
    tot_acc_accr_std_dev = []

  series_infile = open(series_file_names[it], "r")
  try:
    profs_infile = open(profs_file_names[it], "r")
  except:
    print "Could not open the profiles file"

#  print series_file_names[it]
#  print profs_file_names[it]

# calc mean surf precip
  surf_precip = read_my_var(series_infile, "surf_precip")
  mean_surf_precip.append(np.mean(surf_precip[series_from_it:series_to_it]))

# calc acc surf precip
  acc_surf_precip = read_my_var(series_infile, "acc_precip")
  try:
    acc_surf_precip_std_dev = read_my_var(series_infile, "acc_precip_std_dev")
  except:
    print "Could not find acc_precip_std_dev, setting to 0"
    acc_surf_precip_std_dev = np.zeros(len(acc_surf_precip))

  tot_acc_surf_precip.append(acc_surf_precip[series_to_it] - acc_surf_precip[series_from_it])
  tot_acc_surf_precip_std_dev.append(acc_surf_precip_std_dev[series_to_it] + acc_surf_precip_std_dev[series_from_it])

# find cloud base
  try:
    ql = read_my_var(profs_infile, "rliq")
    clbase = np.argmax(ql>qlimit)
  except:
    clbase=0
    print "Could not find the rliq profile, setting clbase=0"
  print clbase

# --- get prflux at cloud base; divide by cloud fraction at this height to get an estimate of average over cloud cells only ---
# TODO: no sense to divide by cloud fraction in a specific cell at a height! We should divide by cloud cover (fraction of cloudy columns), which is in the series file
#  clfrac_at_cbase = read_my_var(profs_infile, "clfrac")[clbase]
#  prflux_at_cbase = read_my_var(profs_infile, "prflux")[clbase]
#  print "cloud fraction at cloud base (calculated using the threshold of LWP in a column from Ackerman et al.): ", clfrac_at_cbase
#  prfluxDivByClFrac.append(read_my_var(profs_infile, "prflux")[clbase] / read_my_var(profs_infile, "clfrac")[clbase])
#  prfluxDivByClFrac_std_dev.append(read_my_var(profs_infile, "prflux_std_dev")[clbase] / read_my_var(profs_infile, "clfrac")[clbase])
# --- get prflux at cloud base height ---
  try:
    prflux.append(read_my_var(profs_infile, "prflux")[clbase] / 2264.705 * 3.6 * 24 ) # convert from W/m2 to mm/day
    prflux_std_dev.append(read_my_var(profs_infile, "prflux_std_dev")[clbase] / 2264.705 * 3.6 * 24 )
  except:
    print "Could not find the prflux profile, setting prflux at cloud base = 0"
    prflux.append(0)
    prflux_std_dev.append(0)
# --- get prflux at cloud base from the prflux vs cloud height profile ---
  try:
    clb_prflux = read_my_var(profs_infile, "base_prflux_vs_clhght")
    clb_prflux_std_dev = np.nan_to_num(read_my_var(profs_infile, "base_prflux_vs_clhght_std_dev"))
    clb_prflux_occur = read_my_var(profs_infile, "base_prflux_vs_clhght number of occurances")
    prfluxFromPrfluxVsClhght.append(np.sum(clb_prflux * clb_prflux_occur) / np.sum(clb_prflux_occur) / 2264.705 * 3.6 * 24 )
    prfluxFromPrfluxVsClhght_std_dev.append(np.sum(clb_prflux_std_dev * clb_prflux_occur) / np.sum(clb_prflux_occur) / 2264.705 * 3.6 * 24 )
  except:
    print "Could not find the base_prflux_vs_clhght data, setting prflux at cloud base from the prflux vs cloud height profile = 0"
    prfluxFromPrfluxVsClhght.append(0)
    prfluxFromPrfluxVsClhght_std_dev.append(0)


# get acnv rate in cloudy cells averaged over whole simulation
  try:
    acc_acnv = read_my_var(series_infile, "acc_cl_acnv25_dycoms")
  except:
    try:
      acc_acnv = read_my_var(series_infile, "acc_cl_acnv25_rico")
    except:
      print "Could not find acc_acnv, setting to 0"
      acc_acnv = np.zeros(len(acc_surf_precip))

  try:
    acc_acnv_std_dev = read_my_var(series_infile, "acc_cl_acnv25_dycoms_std_dev")
  except:
    try:
      acc_acnv_std_dev = read_my_var(series_infile, "acc_cl_acnv25_rico_std_dev")
    except:
      print "Could not find acc_acnv_std_dev, setting to 0"
      acc_acnv_std_dev = np.zeros(len(acc_acnv))

  # acc_cl_acnv25_dycoms_std_dev is a time series of the acnv rate averaged from the start up to the moment.
  # calculate acnv rate averaged from series_from to series_to
  tot_acc_acnv.append(((acc_acnv[series_to_it] * series_to_it) - (acc_acnv[series_from_it] * series_from_it)) / (series_to_it - series_from_it))
  tot_acc_acnv_std_dev.append(acc_acnv_std_dev[series_to_it] + acc_acnv_std_dev[series_from_it])

# get accr rate in cloudy cells averaged over whole simulation
  try:
    acc_accr = read_my_var(series_infile, "acc_cl_accr25_dycoms")
  except:
    try:
      acc_accr = read_my_var(series_infile, "acc_cl_accr25_rico")
    except:
      print "Could not find acc_accr, setting to 0"
      acc_accr = np.zeros(len(acc_surf_precip))

  try:
    acc_accr_std_dev = read_my_var(series_infile, "acc_cl_accr25_dycoms_std_dev")
  except:
    try:
      acc_accr_std_dev = read_my_var(series_infile, "acc_cl_accr25_rico_std_dev")
    except:
      print "Could not find acc_accr_std_dev, setting to 0"
      acc_accr_std_dev = np.zeros(len(acc_accr))

  # acc_cl_accr25_dycoms_std_dev is a time series of the accr rate averaged from the start up to the moment.
  # calculate accr rate averaged from series_from to series_to
  tot_acc_accr.append(((acc_accr[series_to_it] * series_to_it) - (acc_accr[series_from_it] * series_from_it)) / (series_to_it - series_from_it))
  tot_acc_accr_std_dev.append(acc_accr_std_dev[series_to_it] + acc_accr_std_dev[series_from_it])

  if((it+1) % 4 == 0):
   # tot_acc_surf_precip_std_dev = [3 * x for x in tot_acc_surf_precip_std_dev] # we show errors bars with 3 std dev
    tot_acc_surf_precip = [(24. / averaging_period) * x for x in tot_acc_surf_precip] # turn into mm / day
    tot_acc_surf_precip_std_dev = [(24. / averaging_period) * x for x in tot_acc_surf_precip_std_dev] # turn into mm / day
    print tot_acc_surf_precip
    print "prflux at cloud base altitude: ",prflux
    #print "prflux at cloud base altitude divided by cloud fraction: ",prfluxDivByClFrac
    print "prflux at cloud base from the prflux vs cloud height profile: ",prfluxFromPrfluxVsClhght
    #axarr[0].plot(GCCN_conc, mean_surf_precip, 'o')
    tot_acc_acnv = [24. * 3600. * x for x in tot_acc_acnv] # turn into g / m^3 / day
    tot_acc_acnv_std_dev = [24. * 3600. * x for x in tot_acc_acnv_std_dev] # same
    tot_acc_accr = [24. * 3600. * x for x in tot_acc_accr] # turn into g / m^3 / day
    tot_acc_accr_std_dev = [24. * 3600. * x for x in tot_acc_accr_std_dev] # same
    axarr[0,0].errorbar(GCCN_conc, tot_acc_surf_precip, yerr = tot_acc_surf_precip_std_dev, marker='o', fmt='.', label = varlabels[(it)/4])
    axarr[0,1].errorbar(GCCN_conc, prflux, yerr = prflux_std_dev, marker='o', fmt='.')
    axarr[1,0].errorbar(GCCN_conc, tot_acc_acnv, yerr = tot_acc_acnv_std_dev, marker='o', fmt='.')
    axarr[1,1].errorbar(GCCN_conc, tot_acc_accr, yerr = tot_acc_acnv_std_dev, marker='o', fmt='.')

axarr[0,0].set_ylabel('surface precipitation [mm/day]')
axarr[0,1].set_ylabel('cloud base precipitation [mm/day]')
axarr[1,0].set_ylabel('autoconversion [g/(m$^3$ day)]')
axarr[1,1].set_ylabel('accretion [g/(m$^3$ day)]')

#axarr[1,0].set_xlabel('GCCN concentration [cm$^{-3}$]')
#axarr[1,1].set_xlabel('GCCN concentration [cm$^{-3}$]')
axarr[1,0].set_xlabel('$N^\mathrm{init}_\mathrm{GCCN}$ [cm$^{-3}$]')
axarr[1,1].set_xlabel('$N^\mathrm{init}_\mathrm{GCCN}$ [cm$^{-3}$]')

# legend font size
#plt.rcParams.update({'font.size': 10})

# hide axes on empty plots
#if nplots % nploty == 0:
#  nemptyplots=0
#else:
#  nemptyplots = nploty - nplots % nploty
#emptyplots = np.arange(nploty - nemptyplots, nploty)
#for empty in emptyplots:
#  axarr[nplotx-1, empty].axis('off')

#axes = plt.gca()
#axes.tick_params(direction='in')
x_arr = np.arange(nplotx)
y_arr = np.arange(nploty)
for x in x_arr:
  for y in y_arr:
    #tics inside
    axarr[x,y].tick_params(direction='in', which='both', top=1, right=1)
    #minor tics
    axarr[x,y].xaxis.set_minor_locator(AutoMinorLocator())
    axarr[x,y].yaxis.set_minor_locator(AutoMinorLocator())
    #labels and tics font size
    for item in ([axarr[x,y].xaxis.label, axarr[x,y].yaxis.label] + axarr[x,y].get_xticklabels() + axarr[x,y].get_yticklabels()):
      item.set_fontsize(10)
    # subplot numbering
    axarr[x,y].text(0.1, 0.9, labeldict[x*nploty + y], fontsize=10, transform=axarr[x,y].transAxes)

## show legends
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x].legend(loc="upper center")

#single legend for the whole figure
handles, labels = axarr[0,0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.475,0))


#figure size
fig.set_size_inches(5.5, 5.5)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.1 + (len(labels) - 2) * 0.1, wspace=0.3)#, hspace=0.25)

#fig.tight_layout(pad=0.3, w_pad=0, h_pad=0)

#figure size
#fig.set_size_inches(7.874, 6 + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0.15 + (len(labels) - 2) * 0.02, hspace=0.25)


#plt.show()
fig.savefig(sys.argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))

