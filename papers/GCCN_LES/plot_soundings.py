import numpy as np
import matplotlib.pyplot as plt

#plt.rcParams.update({'font.size': 40})
#plt.figure(figsize=(40,40))

plt.figure(figsize=(3.2,2.5))

fig, axarr = plt.subplots(1, 4 )

#stratocumulus
h = 1500
zi = 795

z = np.arange(0,h,1, dtype=float)
thl = z.copy()
qt = z.copy()

thl[z<zi] = 288.3
thl[z>=zi] = 295 + (z[z>=zi] - zi)**(1./3.)

qt[z<zi] = 9.45
print(np.exp(-(z[z>=zi] - zi) / 500.))
qt[z>=zi] = 5 - 3 * (1. - np.exp(-(z[z>=zi] - zi) / 500.))

axarr[0].plot(thl,z)
axarr[1].plot(qt,z)

axarr[0].set_ylabel('height [m]')
axarr[0].set_xlabel(r'$\theta_l$ [K]')
axarr[1].set_xlabel(r'$q_t$ [g/kg]')

axarr[0].set_title('stratocumulus')
axarr[1].set_title('stratocumulus')

#cumulus
h = 4000

z = np.arange(0,h,1, dtype=float)
th = z.copy()
qv = z.copy()

zi = 740.
th[z<zi] = 297.9
th[z>=zi] = 297.9 + (z[z>=zi] - zi) / (4000. - 740.)*(317. - 297.9)


zi = 3260
qv[z>=zi] = 2.4 + (z[z>=zi] - zi) / (4000. - 3260.) *(1.8-2.4)
qv[np.logical_and(z>=740, z<3260)] = 13.8 + (z[np.logical_and(z>=740, z<3260)] - 740) / (3260. - 740.) *(2.4 - 13.8)
zi = 740
qv[z<zi] = 16. + (z[z<zi]) / (740.)*(13.8 - 16.)

axarr[2].plot(th,z)
axarr[3].plot(qv,z)
axarr[2].set_xlabel(r'$\theta$ [K]')
axarr[3].set_xlabel(r'$q_v$ [g/kg]')

axarr[2].set_title('cumulus')
axarr[3].set_title('cumulus')

fig.tight_layout(pad=0.3, w_pad=0, h_pad=0)

axarr[0].set_ylim(0,1500)
axarr[1].set_ylim(0,1500)
axarr[2].set_ylim(0,4000)
axarr[3].set_ylim(0,4000)

axarr[1].set_yticklabels([])
axarr[3].set_yticklabels([])

plt.savefig('/home/piotr/praca/moje_publikacje/GCCN/paper/figs/initial_soundings.pdf')




