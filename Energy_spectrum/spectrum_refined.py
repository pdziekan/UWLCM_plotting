# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

# names of arrays with data to be analyzed (in timestep files)
velocities = ["u", "v", "w"]
velocities = ["w"]

#vel_suffices = ["", " reconstructed", " refined"]
#pos_suffices = ["", " refined", " refined"]

vel_suffices = ["", " reconstructed"]
pos_suffices = ["", " refined"]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
_from_lvl = int(argv[4])
_to_lvl = int(argv[5])

directories = argv[6:len(argv):2]
labels = argv[7:len(argv):2]
print directories, labels

# directories loop
for directory, lab in zip(directories, labels):
  Exy_avg = {}
  #dx = {}
  #ref = {} # refinement = dx / dx_refined
  #suffixes loop, suffix are for resolved/refined/reconstructed data
  for vel_suf, pos_suf in zip(vel_suffices, pos_suffices):
    # read nx, ny, nz, dx, dy, dz
    w3d = h5py.File(directory + "/const.h5", "r")["X"+pos_suf][:,:,:]
    nx, ny, nz = w3d.shape
    dx =  w3d[1][0][0] - w3d[0][0][0]
    # it is asumed that dx == dy
    ref = int(h5py.File(directory + "/const.h5", "r")['advection'].attrs['di'][0] / dx)
    print dx,ref

    # initiliaze average for each velocity
    for _vel in velocities:
      vel = _vel + vel_suf
      Exy_avg[vel] = np.zeros((nx-1)/2 + 1)

#      pos = _pos + pos_suf
#      pos_arr = h5py.File(directory + "/const.h5", "r")[pos]
#      dx[vel] = pos_arr[1][1][1] - pos_arr[0][0][0]
#      print dx[vel]
#      ref[vel_suf] = int(dx[_vel + vel_suffices[0]]/dx[vel])
#      print ref[vel_suf]
#      print fconst.keys()
#      adve_group = fconst["advection"]
#      print adve_group
#      dx[vel] = h5py.File(directory + "/const.h5", "r")[cs]
#      print dx

    from_lvl = _from_lvl * ref
    to_lvl =   _to_lvl * ref
    
    for t in range(time_start, time_end+1, outfreq):
      filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
      print filename
  
      for _vel in velocities:
        vel = _vel + vel_suf
        print vel
    
        print nx,ny
        w3d = h5py.File(filename, "r")[vel][0:nx-1,0:ny-1,:] # * 4. / 3. * 3.1416 * 1e3 
        
        for lvl in range(from_lvl, to_lvl+1):
          w2d = w3d[:, :, lvl]
          #print w2d
          
          wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
          wky = 1.0 / np.sqrt(ny - 1) * np.fft.rfft(w2d, axis = 1)
          
          Ex = 0.5 * (np.abs(wkx) ** 2)
#          Ex = (np.abs(wkx) ** 2)
          Ex = np.mean(Ex, axis = 1)
          Ey = 0.5 * (np.abs(wky) ** 2)
          #Ey = (np.abs(wky) ** 2)
          Ey = np.mean(Ey, axis = 0)
          
#          Exy = 0.5 * (Ex + Ey)
          Exy = Ex
          Exy_avg[vel] += Exy
    
        K = np.fft.rfftfreq(nx - 1) / dx
    #    plt.loglog(K, Exy)
        #lmbd = dx[vel] / K
        lmbd = 1 / K
        print K, lmbd
      
      if (t == time_start and lab==labels[0]):
        plt.loglog(lmbd, 2e-6* K**(-5./3.) )
    
    for _vel in velocities:
      vel = _vel + vel_suf
      Exy_avg[vel] /= (time_end - time_start) / outfreq + 1
      Exy_avg[vel] /= to_lvl+1 - from_lvl
    #crudely scale
    #  Exy_avg[vel] /= Exy_avg[vel][len(Exy_avg[vel])-1]
      plt.loglog(lmbd, Exy_avg[vel] , linewidth=2, label=lab+"_"+vel)
 
#plt.xlim(10**4,10**2)
plt.gca().invert_xaxis()
plt.xlabel("l[m]")
plt.ylabel("E")
plt.legend()
plt.grid(True, which='both', linestyle='--')
#plt.title("Mean PSD of w 322m<z<642m @3h")
plt.show()
