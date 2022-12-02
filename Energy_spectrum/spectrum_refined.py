# (C) Maciej Waruszewski, Piotr Dziekan

import argparse
import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = 10, 10

parser = argparse.ArgumentParser(description='Plot energy spectra from UWLCM output.')

parser.add_argument("-v", "--vars", action="extend", nargs="+", type=str, help="list of variables to be plotted", required=True)
parser.add_argument("-ts", "--time_start", type=int, required=True, help="start of the averaging period [s]")
parser.add_argument("-te", "--time_end", type=int, required=True, help="end of the averaging period [s]")
parser.add_argument("-ls", "--level_start", type=int, required=True, help="lowest level of the averaging area [m]")
parser.add_argument("-le", "--level_end", type=int, required=True, help="highest level of the averaging area [m]")
parser.add_argument("-d", "--dirs", action="extend", nargs="+", type=str, help="list of directories with the data", required=True)
parser.add_argument("-l", "--labels", action="extend", nargs="+", type=str, help="list of labels of the data (same order as --dirs)", required=True)
parser.add_argument("-of", "--outfig", help="output file name", required=True)
args = parser.parse_args()


# names of arrays with data to be analyzed (in timestep files)
#args.vars = ["u", "v", "w"]
#args.vars = ["w"]

#var_suffices = ["", " reconstructed", " refined"]
#pos_suffices = ["", " refined", " refined"]

#var_suffices = ["", " reconstructed"]
#pos_suffices = ["", " refined"]

#TODO: plotting refined arrays still needs work
var_suffices = [""]
pos_suffices = [""]

#time_start = int(argv[1])
#time_end = int(argv[2])
#outfreq = int(argv[3])
#_from_lvl = int(argv[4])
#_to_lvl = int(argv[5])

#directories = argv[6:len(argv):2]
#labels = argv[7:len(argv):2]
#print((directories, labels))

# directories loop
for directory, lab in zip(args.dirs, args.labels):
  Exy_avg = {}

  # read some constant parameters
  with h5py.File(directory + "/const.h5", 'r') as consth5:
    user_params = consth5.get("user_params")
    outfreq = int(user_params.attrs["outfreq"][0])
    advection = consth5.get("advection")
    dx_adve = advection.attrs["di"] # its the resolved dx
    dz_adve = advection.attrs["dk"] # its the resolved dx
    dt = advection.attrs["dt"]

  time_start_idx = int(args.time_start / dt)
  time_end_idx = int(args.time_end / dt)
  level_start_idx = int(args.level_start / dz_adve)
  level_end_idx = int(args.level_end / dz_adve)


  #dx = {}
  #ref = {} # refinement = dx / dx_refined
  #suffixes loop, suffix are for resolved/refined/reconstructed data
  for var_suf, pos_suf in zip(var_suffices, pos_suffices):
    # read nx, ny, nz, dx, dy, dz
    w3d = h5py.File(directory + "/const.h5", "r")["X"+pos_suf][:,:,:]
    nx, ny, nz = tuple(x-1 for x in w3d.shape)
    dx =  w3d[1][0][0] - w3d[0][0][0]
    # it is asumed that dx == dy
    ref = int(dx_adve / dx)
    print((nx,dx,ref))

    # initiliaze average for each velocity
    for _var in args.vars:
      var = _var + var_suf
      Exy_avg[var] = np.zeros(int((nx-1)/2 + 1))

#      pos = _pos + pos_suf
#      pos_arr = h5py.File(directory + "/const.h5", "r")[pos]
#      dx[var] = pos_arr[1][1][1] - pos_arr[0][0][0]
#      print dx[var]
#      ref[var_suf] = int(dx[_var + var_suffices[0]]/dx[var])
#      print ref[var_suf]
#      print fconst.keys()
#      adve_group = fconst["advection"]
#      print adve_group
#      dx[var] = h5py.File(directory + "/const.h5", "r")[cs]
#      print dx

    from_lvl = level_start_idx * ref
    to_lvl =   level_end_idx * ref
    
    for t in range(time_start_idx, time_end_idx+1, outfreq):
      filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
      print(filename)
  
      for _var in args.vars:
        var = _var + var_suf
        print(var)
    
        print(nx,ny)
        w3d = h5py.File(filename, "r")[var][0:nx-1,0:ny-1,:] # * 4. / 3. * 3.1416 * 1e3 
        
        for lvl in range(from_lvl, to_lvl+1):
          w2d = w3d[:, :, lvl]
          #print w2d
          
          wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
          wky = 1.0 / np.sqrt(ny - 1) * np.fft.rfft(w2d, axis = 1)
          
#          Ex = 0.5 * (np.abs(wkx) ** 2)
          Ex = (np.abs(wkx) ** 2)
          Ex = np.mean(Ex, axis = 1)
#          Ey = 0.5 * (np.abs(wky) ** 2)
          Ey = (np.abs(wky) ** 2)
          Ey = np.mean(Ey, axis = 0)
          
          Exy = 0.5 * (Ex + Ey)
#          Exy = Ex
          Exy_avg[var] += Exy
    
        K = np.fft.rfftfreq(nx - 1) / dx
    #    plt.loglog(K, Exy)
        #lmbd = dx[var] / K
        lmbd = 1 / K
        print(K, lmbd)
      
      if (t == time_start_idx and lab==args.labels[0]):
        plt.loglog(lmbd, 2e-6* K**(-5./3.) )
    
    for _var in args.vars:
      var = _var + var_suf
      Exy_avg[var] /= (time_end_idx - time_start_idx) / outfreq + 1
      Exy_avg[var] /= to_lvl+1 - from_lvl
      #crudely scale
      #Exy_avg[var] /= Exy_avg[var][len(Exy_avg[var])-1]
      Exy_avg[var] /= np.sum(Exy_avg[var])
      plt.loglog(lmbd, Exy_avg[var] , linewidth=2, label=lab+"_"+var)
 
#plt.xlim(10**4,10**2)
plt.gca().invert_xaxis()
plt.xlabel("l[m]")
plt.ylabel("E")
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.title("Mean PSD z=["+str(args.level_start)+"m, "+str(args.level_end)+"m] @["+str(args.time_start)+"s, "+str(args.time_end)+"s]")
plt.savefig(args.outfig)
plt.show()
