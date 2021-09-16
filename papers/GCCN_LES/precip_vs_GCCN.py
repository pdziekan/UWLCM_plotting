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
relative = int(sys.argv[6])
cloud_type = int(sys.argv[7]) # ==0 - stratocumulus, !=0 - cumulus

#varlabels = ["{\it ScNc30}", "{\it ScNc40\_salt\_CCN}", "{\it ScNc45}", "{\it ScNc105}"]

if cloud_type != 0:
  varlabels = ["{\it Cu38}", "{\it Cu60}", "{\it Cu88}"]
  CCN_conc = [52.5, 105, 210]
  nc = [38, 60, 85]
  ensemble = [5,6,4]
  cl_cover_series_name = "cloud_cover_rico"
else:
  varlabels = ["{\it Sc38}", "{\it Sc59}", "{\it Sc115}"]
  CCN_conc = [95, 190, 475]
  nc = [38, 59, 115]
  ensemble = [4,2,2]
  cl_cover_series_name = "cloud_cover_dycoms"

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
file_no = np.arange(8, len(sys.argv)-1 , 1)
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
    tot_cl_cover_mean = []
    tot_cl_cover_std_dev = []
    mean_surf_precip = []
    tot_acc_surf_precip = []
    tot_acc_surf_precip_std_dev = []
    tot_acc_surf_precip_over_cl_cover = []
    tot_acc_surf_precip_over_cl_cover_std_dev = []
    prflux = []
    prflux_std_dev = []
    prflux_over_cl_cover = []
    prflux_over_cl_cover_std_dev = []
    tot_prflux_over_tot_cl_cover_mean = []
    tot_prflux_over_tot_cl_cover_std_dev = []
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
  surf_precip_std_dev = read_my_var(series_infile, "surf_precip_std_dev")
  mean_surf_precip.append(np.mean(surf_precip[series_from_it:series_to_it]))

# calc acc surf precip
  acc_surf_precip = read_my_var(series_infile, "acc_precip")
  try:
    acc_surf_precip_std_dev = read_my_var(series_infile, "acc_precip_std_dev")
  except:
    print "Could not find acc_precip_std_dev, setting to 0"
    acc_surf_precip_std_dev = np.zeros(len(acc_surf_precip))

# read cloud cover
  cl_cover = read_my_var(series_infile, cl_cover_series_name)
  try:
    cl_cover_std_dev = read_my_var(series_infile, cl_cover_series_name+"_std_dev")
  except:
    print "Could not find acc_precip_std_dev, setting to 0"
    cl_cover_std_dev = np.zeros(len(cl_cover))

  tot_cl_cover_mean.append(np.mean(cl_cover[series_from_it:series_to_it]))
  tot_cl_cover_std_dev.append(np.sqrt(np.sum(np.power(cl_cover_std_dev[series_from_it:series_to_it],2)))) # propagation of uncertainty (wikipedia)

  tot_acc_surf_precip.append(acc_surf_precip[series_to_it] - acc_surf_precip[series_from_it])
  tot_acc_surf_precip_std_dev.append(acc_surf_precip_std_dev[series_to_it] + acc_surf_precip_std_dev[series_from_it])
 
  surf_precip_over_cl_cover = surf_precip / cl_cover
  surf_precip_over_cl_cover_std_dev = surf_precip / cl_cover * np.sqrt(np.power(surf_precip_std_dev / surf_precip,2) + np.power(cl_cover_std_dev / cl_cover,2))
  tot_acc_surf_precip_over_cl_cover.append(np.mean(surf_precip_over_cl_cover[series_from_it:series_to_it]))
  tot_acc_surf_precip_over_cl_cover_std_dev.append(np.mean(surf_precip_over_cl_cover_std_dev[series_from_it:series_to_it]))
#  tot_acc_surf_precip_over_cl_cover_std_dev = []

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
# --- prflux at cloud base height divided by average cloud cover---
  prflux_over_cl_cover.append(prflux[-1] / np.mean(cl_cover[series_from_it:series_to_it]))
  prflux_over_cl_cover_std_dev.append(prflux_over_cl_cover[-1] * np.sqrt(np.power(prflux_std_dev[-1] / prflux[-1],2) + np.power(np.mean(cl_cover_std_dev[series_from_it:series_to_it]) / np.mean(cl_cover[series_from_it:series_to_it]),2)))

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
    tot_acc_acnv = [24. * 3600. * x for x in tot_acc_acnv] # turn into g / m^3 / day
    tot_acc_acnv_std_dev = [24. * 3600. * x for x in tot_acc_acnv_std_dev] # same
    tot_acc_accr = [24. * 3600. * x for x in tot_acc_accr] # turn into g / m^3 / day
    tot_acc_accr_std_dev = [24. * 3600. * x for x in tot_acc_accr_std_dev] # same

    tot_acc_surf_precip_over_tot_cl_cover_mean = [x / y for x,y in zip(tot_acc_surf_precip, tot_cl_cover_mean)]
    tot_acc_surf_precip_over_tot_cl_cover_std_dev = [x / y * np.sqrt(np.power(sx/x,2) + np.power(sy/y,2)) for x,sx,y,sy in zip(tot_acc_surf_precip, tot_acc_surf_precip_std_dev, tot_cl_cover_mean, tot_cl_cover_std_dev)]

    tot_prflux_over_tot_cl_cover_mean = [x / y for x,y in zip(prflux, tot_cl_cover_mean)]
    tot_prflux_over_tot_cl_cover_std_dev = [x / y * np.sqrt(np.power(sx/x,2) + np.power(sy/y,2)) for x,sx,y,sy in zip(prflux, prflux_std_dev, tot_cl_cover_mean, tot_cl_cover_std_dev)]
    
#    tot_acc_surf_precip / tot_cl_cover_mean * np.sqrt(np.power(tot_acc_surf_precip_std_dev / tot_acc_surf_precip,2) + np.power(tot_cl_cover_std_dev / tot_cl_cover_mean,2))

    print "surf precip", tot_acc_surf_precip
    print "surf precip divided by cloud cover", tot_acc_surf_precip_over_cl_cover
    print "prflux at cloud base altitude: ",prflux
    #print "prflux at cloud base altitude divided by cloud fraction: ",prfluxDivByClFrac
    print "prflux at cloud base from the prflux vs cloud height profile: ",prfluxFromPrfluxVsClhght
    print "prflux at cloud base altitude divided by cloud cover: ",prflux_over_cl_cover
    print "acnv", tot_acc_acnv
    print "accr", tot_acc_accr
    print "tot_acc_surf_precip_over_tot_cl_cover_mean", tot_acc_surf_precip_over_tot_cl_cover_mean
    print "tot_acc_surf_precip_over_tot_cl_cover_std_dev", tot_acc_surf_precip_over_tot_cl_cover_std_dev
    #axarr[0].plot(GCCN_conc, mean_surf_precip, 'o')
    if(relative):
      GCCN_CCN_rat = [GCCN_con / CCN_conc[int(np.floor(it/4))] for GCCN_con in GCCN_conc]
      #GCCN_CCN_rat = [GCCN_con / nc[int(np.floor(it/4))] for GCCN_con in GCCN_conc]
      axarr[0,0].errorbar(GCCN_CCN_rat, tot_acc_surf_precip - tot_acc_surf_precip[0], yerr = (tot_acc_surf_precip_std_dev + tot_acc_surf_precip_std_dev[0]) / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.', label = varlabels[(it)/4])
      axarr[0,1].errorbar(GCCN_CCN_rat, prflux - prflux[0]                          , yerr = (prflux_std_dev + prflux_std_dev[0])                           / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
      axarr[1,0].errorbar(GCCN_CCN_rat, tot_acc_acnv - tot_acc_acnv[0]              , yerr = (tot_acc_acnv_std_dev + tot_acc_acnv_std_dev[0])               / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
      axarr[1,1].errorbar(GCCN_CCN_rat, tot_acc_accr - tot_acc_accr[0]              , yerr = (tot_acc_acnv_std_dev + tot_acc_acnv_std_dev[0])               / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
    else:
#      axarr[0,0].errorbar(GCCN_conc, tot_acc_surf_precip,               yerr = tot_acc_surf_precip_std_dev / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.', label = varlabels[(it)/4])
      axarr[0,0].errorbar(GCCN_conc, tot_acc_surf_precip_over_tot_cl_cover_mean,  yerr = tot_acc_surf_precip_over_tot_cl_cover_std_dev / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.', label = varlabels[(it)/4])
#      axarr[0,0].errorbar(GCCN_conc, tot_acc_surf_precip_over_cl_cover, yerr = tot_acc_surf_precip_over_cl_cover_std_dev / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.', label = varlabels[(it)/4])
#      axarr[0,1].errorbar(GCCN_conc, prflux,                            yerr = prflux_std_dev              / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
      axarr[0,1].errorbar(GCCN_conc, tot_prflux_over_tot_cl_cover_mean,           yerr = tot_prflux_over_tot_cl_cover_std_dev / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
#      axarr[0,1].errorbar(GCCN_conc, prflux_over_cl_cover,              yerr = prflux_over_cl_cover_std_dev / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
      axarr[1,0].errorbar(GCCN_conc, tot_acc_acnv,                      yerr = tot_acc_acnv_std_dev        / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')
      axarr[1,1].errorbar(GCCN_conc, tot_acc_accr,                      yerr = tot_acc_acnv_std_dev        / np.sqrt(ensemble[(it)/4]), marker='o', fmt='.')

if(relative):
  axarr[0,0].set_ylabel('exc. surf. precip. [mm/day]')
  axarr[0,1].set_ylabel('exc. cl. base precip. [mm/day]')
  axarr[1,0].set_ylabel('exc. autoconv. [g/(m$^3$ day)]')
  axarr[1,1].set_ylabel('exc. accr. [g/(m$^3$ day)]')
  #axarr[1,0].set_xlabel('GCCN concentration [cm$^{-3}$]')
  #axarr[1,1].set_xlabel('GCCN concentration [cm$^{-3}$]')
  axarr[1,0].set_xlabel('$N^\mathrm{init}_\mathrm{GCCN} / N_\mathrm{CCN}$')
  axarr[1,1].set_xlabel('$N^\mathrm{init}_\mathrm{GCCN} / N_\mathrm{CCN}$')
else:
  axarr[0,0].set_ylabel('$C_\mathrm{surf}$ [mm/day]')
  axarr[0,1].set_ylabel('$C_\mathrm{clb}$ [mm/day]')
  axarr[1,0].set_ylabel('$R_\mathrm{acnv}$ [g/(m$^3$ day)]')
  axarr[1,1].set_ylabel('$R_\mathrm{accr}$ [g/(m$^3$ day)]')
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

