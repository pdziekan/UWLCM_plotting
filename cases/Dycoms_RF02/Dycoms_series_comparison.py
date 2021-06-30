from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")

from Dycoms_reference_plots import plot_reference_series
from plot_ranges import xscaledict, yscaledict, xlimdict_series, ylimdict_series
from plot_series import *
from latex_labels import labeldict
from autoscale_y import *

# activate latex text rendering
rc('text', usetex=True)

dycoms_vars = ["lwp", "er", "wvarmax", "surf_precip", "cl_nc", "cl_nr", "cloud_base_dycoms", "cloud_cover_dycoms", "cl_acnv25_dycoms", "cl_accr25_dycoms"]# "cloud_base"]#, "cfrac"]
#dycoms_vars = ["clb_bigrain_mean_rd","clb_bigrain_mean_kappa","clb_bigrain_mean_conc","clb_bigrain_mean_inclt", "cl_nr"]

# variables that need rescaling of the yrange to the limited x range of 1-6h
rescale_vars = ["lwp", "er", "wvarmax", "cl_nc", "cloud_base_dycoms", "cloud_cover_dycoms"]

# init the plot
nplotx = 3
nploty= 4
x_arr = np.arange(nplotx)
y_arr = np.arange(nploty)
fig, axarr = plt.subplots(nplotx,nploty)

for x in x_arr:
  for y in y_arr:
    axarr[x,y].margins(y=0.2)

#plot_reference_series(dycoms_vars, 0, nplotx, nploty, axarr)
plot_series(dycoms_vars, 0, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_series, ylimdict_series, xlabel='Time [h]')

# show legends on each subplot
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x,y].legend()

# legend font size
plt.rcParams.update({'font.size': 8})

# hide axes on empty plots
if len(dycoms_vars) % nploty == 0:
  nemptyplots = 0
else:
  nemptyplots = nploty - len(dycoms_vars) % nploty
emptyplots = np.arange(nploty - nemptyplots, nploty)
for empty in emptyplots:
  axarr[nplotx-1, empty].axis('off')

#axes = plt.gca()
#axes.tick_params(direction='in')
for x in x_arr:
  for y in y_arr:
    #tics inside
    axarr[x,y].tick_params(direction='in', which='both', top=1, right=1)
    #minor tics
    axarr[x,y].xaxis.set_minor_locator(AutoMinorLocator())
    axarr[x,y].yaxis.set_minor_locator(AutoMinorLocator())
    #labels and tics font size
    for item in ([axarr[x,y].xaxis.label, axarr[x,y].yaxis.label] + axarr[x,y].get_xticklabels() + axarr[x,y].get_yticklabels()):
      item.set_fontsize(8)
    # subplot numbering
    if y < nploty - nemptyplots or x < (nplotx - 1): #nonempty plots
      axarr[x,y].text(0.2, 0.875, labeldict[y + x*nploty], fontsize=8, transform=axarr[x,y].transAxes)

      # rescale y range to the visible x range, note: overrides ylim!
      var = dycoms_vars[x*nploty + y]
      if var in rescale_vars:
        autoscale_y(axarr[x,y], margin=0.3)
      
      # hide hrzntl tic labels
      if x*nploty + y < nplotx * nploty - nemptyplots - nploty:
        axarr[x,y].set_xticklabels([])

#single legend for the whole figure
handles, labels = axarr[0,0].get_legend_handles_labels()

lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))

#figure size
fig.set_size_inches(7.874, 5. + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.18 + (len(labels) - 2) * 0.03, hspace=0.1, wspace=0.4)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)

#plt.show()
fig.savefig(argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))
