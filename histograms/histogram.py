
import argparse
import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import OrderedDict

c_pd = 1005.7
R_d = 287.

def exner(p):
  return (p / 1e5)**(R_d/c_pd)
v_exner = np.vectorize(exner)

def calc_T(th, p):
  return th * v_exner(p)
v_calc_T = np.vectorize(calc_T)

# Tetens: r_vs=3.8/(p*exp(-17.2693882*(T-273.15)/(T-35.86))-6.109)  p in mb, T in Kelvins
def calc_rv_s(th, rv, p):
  T = v_calc_T(th, p)
  return 3.8 / (p*np.exp(-17.2693882*(T-273.15)/(T-35.86))-6.109)
v_calc_rv_s = np.vectorize(calc_rv_s)

def calc_RH(th, rv, p):
  return rv / v_calc_rv_s(th, rv, p) / 100.
v_calc_RH = np.vectorize(calc_RH)

#mpl.rcParams['figure.figsize'] = 10, 10
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(10,10))

parser = argparse.ArgumentParser(description='Plot histograms of variables from UWLCM output.')

parser.add_argument("-v", "--vars", action="extend", nargs="+", type=str, help="list of variables to be plotted", required=True)
parser.add_argument("-ts", "--time_start", type=float, required=True, help="start of the averaging period [s]")
parser.add_argument("-te", "--time_end", type=float, required=True, help="end of the averaging period [s]")
parser.add_argument("-ls", "--level_start", type=float, required=True, help="lowest level of the averaging area [m]")
parser.add_argument("-le", "--level_end", type=float, required=True, help="highest level of the averaging area [m]")
parser.add_argument("-d", "--dirs", action="extend", nargs="+", type=str, help="list of directories with the data", required=True)
parser.add_argument("-l", "--labels", action="extend", nargs="+", type=str, help="list of labels of the data (same order as --dirs)", required=True)
parser.add_argument("-of", "--outfig", help="output file name", required=True)
parser.add_argument("--outfreq", type=int, required=False, help="output frequency of the simulation [number of time steps], if not specified it will be read from const.h5 (if possible)")
parser.add_argument('--normalize', action='store_true', help="normalize the histogram")
parser.add_argument('--mask_rico', action='store_true', help="compute histogram only within cloud cells (using the rico cloud mask)")
parser.set_defaults(normalie=False)
args = parser.parse_args()


nx = {}
ny = {}
nz = {}
dx = {}
dz = {}
ref = {}


# directories loop
for directory, lab in zip(args.dirs, args.labels):
  print(directory, lab)
  can_plot_refined_RH_derived = True

  total_arr   = OrderedDict()
  plot_labels = OrderedDict()

  # init parameters from const.h5
  with h5py.File(directory + "/const.h5", 'r') as consth5:
    user_params = consth5.get("user_params")
    if args.outfreq is None:
      outfreq = int(user_params.attrs["outfreq"][0])
    else:
      outfreq = args.outfreq
    advection = consth5.get("advection")
    dx_adve = advection.attrs["di"] # its the resolved dx
    dz_adve = advection.attrs["dk"] # its the resolved dx
    dt = advection.attrs["dt"]
    nx_adve = consth5["X"][:,:,:].shape[0] - 1
    nz_adve = consth5["Z"][:,:,:].shape[2] - 1
    X = dx_adve * (nx_adve-1)
    Z = dz_adve * (nz_adve-1)
    p_e = consth5["p_e"][:]
    try:
      refined_p_e = consth5["refined p_e"][:]
    except:
      can_plot_refined_RH_derived = False
      print("'refined p_e' not found in const.h5. Won't be able to plot refined_RH_derived")

  # vars loop
  for var in args.vars:
    print(var)

    if(not can_plot_refined_RH_derived and var == "refined RH_derived"):
      continue

    time_start_idx = int(args.time_start / dt)
    time_end_idx = int(args.time_end / dt)


    # init variable-specific array parameters based on the first timestep
    filename = directory + "/timestep" + str(time_start_idx).zfill(10) + ".h5"

    # special case of RH calculated from th and rv
    # NOTE: RH_derived is calculated with r_v / r_vs, where r_vs comes from the Tetens formula
    #       in UWLCM r_vs comes from clausius-clapeyron.
    #       RH_derived is shifted to the right with respect to RH from UWLCM - maybe due to this difference?
    if(var == "RH_derived"):
      w3d = h5py.File(filename, "r")["th"][:,:,:]
    elif(can_plot_refined_RH_derived and var == "refined RH_derived"):
      w3d = h5py.File(filename, "r")["refined th"][:,:,:]
    else:
      w3d = h5py.File(filename, "r")[var][:,:,:]


    nx, ny, nz = tuple(x for x in w3d.shape)
    dx = X / (nx - 1)
    ref = int(dx_adve / dx)
    dz = Z / (nz - 1) 

    print("nx_adve: ", nx_adve)
    print("nx: ", nx)
    print("dx_adve: ", dx_adve)
    print("dx: ", dx)

    print("nz_adve: ", nz_adve)
    print("nz: ", nz)
    print("dz_adve: ", dz_adve)
    print("dz: ", dz)

    print("ref: ", ref)
    assert(float(args.level_start / dz).is_integer())
    assert(float(args.level_end / dz).is_integer())
    level_start_idx = int(args.level_start / dz)
    level_end_idx = int(args.level_end / dz) + 1
    print("level start index for this var: ", level_start_idx)
    print("level end index for this var: ", level_end_idx)

    lab = lab + '_' + str(var)

    total_arr[lab] = np.zeros(0) 
    plot_labels[lab] = lab

    # time loop
    for t in range(time_start_idx, time_end_idx+1, outfreq):
      filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
      print(filename)

      if(var == "RH_derived"):
        th = h5py.File(filename, "r")["th"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        rv = h5py.File(filename, "r")["rv"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        p  = np.empty(th.shape) # 3D pressure array, filled from the 1D profile
        for i in np.arange(p.shape[0]):
          for j in np.arange(p.shape[1]):
            p[i,j] = p_e[level_start_idx:level_end_idx]
        w3d = v_calc_RH(th, rv, p)
        total_arr[lab] = np.append(total_arr[lab], w3d)

      elif(var == "refined RH_derived"):
        th = h5py.File(filename, "r")["refined th"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        rv = h5py.File(filename, "r")["refined rv"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        p  = np.empty(th.shape)
        for i in np.arange(p.shape[0]):
          for j in np.arange(p.shape[1]):
            p[i,j] = refined_p_e[level_start_idx:level_end_idx]
        w3d = v_calc_RH(th, rv, p)
        total_arr[lab] = np.append(total_arr[lab], w3d)

      else:
        w3d = h5py.File(filename, "r")[var][0:nx-1, 0:ny-1, level_start_idx:level_end_idx] # * 4. / 3. * 3.1416 * 1e3 
        total_arr[lab] = np.append(total_arr[lab], w3d)


# convert to typical units
#if data == "rain_rw_mom3":
#  total_arr[data][lab] *= 4./3. * 3.1416 * 1e3 * 1e3 # [g/kg]
#if data == "precip_rate":
#  total_arr[data][lab] *= 4./3. * 3.1416 * 1e3 / CellVol * L_evap

#for lab in labels:
##  print  np.average(total_arr[lab])
#  plot_labels[lab] = plot_labels[lab] + '\n <q_r> = {:.3e}'.format(np.average(total_arr["rain_rw_mom3"][lab])) \
#                                      + '\n <precip flux> = {:.3e}'.format(np.average(total_arr["precip_rate"][lab])) \
#                                      + '\n <cloud base lvl> = {:.2f}'.format(np.average(tot_cloud_base_lvl[lab] * 5))
  #_ = plt.hist(total_arr["rain_rw_mom3"].values(), bins='auto', label=plot_labels.values(), density=True)
  print(total_arr)
  data = list(total_arr.values())
  print(data)
  
  # for lin plots:
  n, bins, patches = plt.hist(data, bins=100, label=plot_labels.values(), density=args.normalize, histtype='step', linewidth=2)
#  plt.axvline(x = np.average(data), ls='--', color=patches[0].get_facecolor()[0:3])
  print("total number of cells = " + str(np.sum(n)))

  # for log plots:
  #_ = plt.hist(data, bins=np.logspace(np.log10(1e-6), np.log10(np.amax(data)), 100), label=plot_labels.values(), density=False, histtype='step', linewidth=6)
 # plt.xscale('log')
 # plt.yscale('log')

plt.legend()#loc = 'lower center')
#plt.xlabel('q_r [g/kg]')
ylabel =  'PDF' if args.normalize else '# of cells'
plt.ylabel(ylabel)
plt.grid(True, which='both', linestyle='--')
plt.title("z=["+str(args.level_start)+"m, "+str(args.level_end)+"m] @["+str(args.time_start)+"s, "+str(args.time_end)+"s]")

#plt.savefig('rain_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
plt.savefig(args.outfig)
plt.show()
