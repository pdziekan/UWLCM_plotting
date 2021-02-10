from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.ticker import FormatStrFormatter, NullFormatter
from matplotlib.ticker import MaxNLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
#sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Dycoms_RF02/")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../cases/RICO11/")

from plot_ranges import * 
from plot_profs import *

# activate latex text rendering
rc('text', usetex=True)
# font size
plt.rcParams.update({'font.size': 10})
# fig size
plt.figure(figsize=(7.2,4))

linestyles = ['-', '-.', '--']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

plt.xlabel('mean wet radius of droplets formed on aerosols with $r_\mathrm{dry}>2\, \mu\mathrm{m}$ [$\mu\mathrm{m}$]')
plt.ylabel('height / inversion height')

for no in [2,4,6]:
  file_name = sys.argv[no] + "profiles_7200_18000.dat"
  profs_file = open(file_name, "r")
  my_hgt = read_my_var(profs_file, "position")
  my_gccn_rw = read_my_var(profs_file, "gccn_rw_cl_down")
  
  profs_file.close()
  
  plt.plot(my_gccn_rw, my_hgt, label=sys.argv[no+1], ls=linestyles[int(no/2-1)])

plt.legend()
plt.ylim(0,1.1)
plt.savefig(sys.argv[1])

plt.clf()
for no in [2,4,6]:
  file_name = sys.argv[no] + "profiles_7200_18000.dat"
  profs_file = open(file_name, "r")
  my_hgt = read_my_var(profs_file, "position")
  my_gccn_rw = read_my_var(profs_file, "gccn_rw_cl_up")
  
  profs_file.close()
  
  plt.plot(my_gccn_rw, my_hgt, label=sys.argv[no+1])

plt.legend()
plt.ylim(0,1.1)
plt.savefig(sys.argv[1]+'updraft.pdf')
