
import argparse
import h5py
import numpy as np
from sys import path
import matplotlib.pyplot as plt
from collections import OrderedDict

path.append('/home/piotr/singu_built_libraries/usr/lib/python3/dist-packages/')
from libcloudphxx import common


#c_pd = 1005.7
#R_d = 287.
#
#def exner(p):
#  return (p / 1e5)**(R_d/c_pd)
#v_exner = np.vectorize(exner)
#
def calc_T(th, p):
  return np.float64(th) * common.exner(np.float64(p))
#v_calc_T = np.vectorize(calc_T)

#
## Tetens: r_vs=3.8/(p*exp(-17.2693882*(T-273.15)/(T-35.86))-6.109)  p in mb, T in Kelvins
#def calc_rv_s(th, rv, p):
#  T = v_calc_T(th, p)
#  return 3.8 / (p*np.exp(-17.2693882*(T-273.15)/(T-35.86))-6.109)

def calc_r_vs_T(T, p):
  return common.r_vs(np.float64(T), np.float64(p))
v_calc_r_vs_T = np.vectorize(calc_r_vs_T)

def calc_r_vs(th, p):
  T = calc_T(th, p)
  return common.r_vs(T, p)
v_calc_r_vs = np.vectorize(calc_r_vs)

def calc_RH(th, rv, p):
  return rv / v_calc_r_vs(th, p)
v_calc_RH = np.vectorize(calc_RH)

#mpl.rcParams['figure.figsize'] = 10, 10
plt.rcParams.update({'font.size': 10})
plt.figure(1, figsize=(10,10)) # histogram plot
plt.figure(2, figsize=(10,10)) # scatter plot

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
parser.add_argument('--no_histogram', action='store_true', help="dont save the histogram plot")
parser.add_argument('--no_scatter', action='store_true', help="dont calculate correlations and dot save the scatter plot")
parser.add_argument('--mask_rico', action='store_true', help="compute histogram only within cloud cells (using the rico cloud mask)")
parser.add_argument('--scatter_saturation', action='store_true', help="plot saturation line on th vs rv scatter plots")
parser.set_defaults(normalize=False, scatter_saturation=False)
args = parser.parse_args()


nx = {}
ny = {}
nz = {}
dx = {}
dz = {}
ref = {}

scatter_saturation_plotted = False


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
      print("Skipping the refined_RH_derived plot")
      continue

    time_start_idx = int(args.time_start / dt)
    time_end_idx = int(args.time_end / dt)


    # init variable-specific array parameters based on the first timestep
    filename = directory + "/timestep" + str(time_start_idx).zfill(10) + ".h5"

    # special case of RH calculated from th and rv
    if(var == "RH_derived"):
      w3d = h5py.File(filename, "r")["th"][:,:,:]
    elif(can_plot_refined_RH_derived and var == "refined RH_derived"):
      w3d = h5py.File(filename, "r")["refined th"][:,:,:]
    else:
      try:
        w3d = h5py.File(filename, "r")[var][:,:,:]
      except:
        continue


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

    lab_var = lab + '_' + str(var)

    if(args.mask_rico):
      try:
        mask = h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:] 
        if(mask.shape != w3d.shape):
          print("Cloud mask shape is different than "+var+" shape. Skipping the plot.")
          continue
      except:
        print("Can't find cloud_rw_mom3 data, so can't use RICO cloud mask. Skipping the plot.")
        continue


    total_arr[lab_var] = np.zeros(0) 
    plot_labels[lab_var] = lab_var

    # time loop
    for t in range(time_start_idx, time_end_idx+1, outfreq):
      filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
      print(filename)

      # read the variable
      if(var == "RH_derived"):
        th = h5py.File(filename, "r")["th"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        rv = h5py.File(filename, "r")["rv"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        p  = np.empty(th.shape) # 3D pressure array, filled from the 1D profile
        for i in np.arange(p.shape[0]):
          for j in np.arange(p.shape[1]):
            p[i,j] = p_e[level_start_idx:level_end_idx]
        w3d = v_calc_RH(th, rv, p)

      elif(var == "refined RH_derived"):
        th = h5py.File(filename, "r")["refined th"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        rv = h5py.File(filename, "r")["refined rv"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx]
        p  = np.empty(th.shape)
        for i in np.arange(p.shape[0]):
          for j in np.arange(p.shape[1]):
            p[i,j] = refined_p_e[level_start_idx:level_end_idx]
        w3d = v_calc_RH(th, rv, p)

      else:
        w3d = h5py.File(filename, "r")[var][0:nx-1, 0:ny-1, level_start_idx:level_end_idx] # * 4. / 3. * 3.1416 * 1e3 

      # read and apply cloud mask
      if(args.mask_rico):
        mask = h5py.File(filename, "r")["cloud_rw_mom3"][0:nx-1, 0:ny-1, level_start_idx:level_end_idx] * 4./3. * 3.1416 * 1e3
        mask = np.where(mask > 1.e-5, 1., 0.)
        w3d = w3d[(mask == 1)]


      total_arr[lab_var] = np.append(total_arr[lab_var], w3d)


  # correlations oefficients ad scatter plots
  print(lab, " correlation coefficients:")
  for var1 in args.vars:
    lab_var1 = lab + '_' + str(var1)
    if(not lab_var1 in total_arr):
      continue
    for var2 in args.vars[args.vars.index(var1)+1:]:
      if(var1 == var2):
        continue
      lab_var2 = lab + '_' + str(var2)
      if(not lab_var2 in total_arr):
        continue
      
      # correlation coefficient
      corr = np.corrcoef(total_arr[lab_var1].flatten(), total_arr[lab_var2].flatten(), rowvar=False)
      print(lab_var1, " " + lab_var2 + " : ", corr)

      # plot
      plt.figure(2)
      plt.scatter(total_arr[lab_var1].flatten(), total_arr[lab_var2].flatten(), marker='.', s=1, label=lab)
      plt.xlabel(var1)
      plt.ylabel(var2)

      # add saturation line on a rv vs th plot
      # TODO: add this line also to refined rv vs refined th plots
      if(args.scatter_saturation and var1 == "th" and var2 == "rv" and args.level_start == args.level_end and not scatter_saturation_plotted):
        print("plotting a saturation line on the scatter plot")
        press = p_e[level_start_idx]
        min_th = total_arr[lab_var1].min()
        max_th = total_arr[lab_var1].max()
        min_T = calc_T(min_th, press)
        max_T = calc_T(max_th, press)
        v_th = np.linspace(min_th, max_th, 100)
        v_T = np.linspace(min_T, max_T, 100)
        v_r_vs = v_calc_r_vs_T(v_T, press)
        plt.plot(v_th, v_r_vs, c='black', ls='--', label='r_vs')
        scatter_saturation_plotted = True



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
  print("data:", data)

  # return to the histogram fig
  plt.figure(1)
  
  # for lin plots:
  n, bins, patches = plt.hist(data, bins=100, label=plot_labels.values(), density=args.normalize, histtype='step', linewidth=2)
  for i in np.arange(len(data)):
    if(len(data) > 1):
      plt.axvline(x = np.average(data[i]), ls='--', color=patches[i][0].get_facecolor()[0:3])
    else:
      plt.axvline(x = np.average(data[i]), ls='--', color=patches[0].get_facecolor()[0:3])
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

if(not args.no_histogram):
  #plt.savefig('rain_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')
  plt.savefig(args.outfig)

# scatter plot
plt.figure(2)
plt.legend()#loc = 'lower center')

plt.grid(True, which='both', linestyle='--')
plt.title("z=["+str(args.level_start)+"m, "+str(args.level_end)+"m] @["+str(args.time_start)+"s, "+str(args.time_end)+"s]")

if(not args.no_scatter):
  plt.savefig(args.outfig + 'scatter.svg')
