###Times
from sys import argv, path, maxsize
import cProfile, pstats
import h5py
import functools
import numpy as np
from itertools import product
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing
import concurrent.futures
import threading
import concurrent.futures
from matplotlib.ticker import ScalarFormatter
import timeit
import pandas as pd
from matplotlib.gridspec import SubplotSpec
plt.rcParams.update({'font.size': 21})

mean_value = np.load('PATH TO mean_value.npy',allow_pickle=True)
st_dev = np.load('PATH TO st_dev.npy',allow_pickle=True)
st_dev_error = np.load('PATH TO st_dev_error.npy',allow_pickle=True)
mean_error = np.load('PATH TO mean_error.npy',allow_pickle=True)

def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
    "Sign sets of subplots with title"
    row = fig.add_subplot(grid)
    row.set_title(f'{title}\n', fontweight='semibold')
    row.set_frame_on(False)
    row.axis('off')

sub_text = {0:'(a)', 1:'(b)', 2:'(c)', 3:'(d)', 4:'(e)', 5:'(f)', 6:'(g)', 7:'(h)'}
timesteps = [0.05, 0.1, 0.5]
subplots_text = ['D','LR','MR','HR']

P = int(len(subplots_text))
fig, ax = plt.subplots(2, P, sharex=True)
fig.set_size_inches(19.5, 16.0)
xlabels = [0.05, 0.1, 0.5]

for j in range(len(subplots_text)):        
    ax[0,j].errorbar(timesteps, np.array_split(mean_value, P)[j],  color='k',  yerr=np.array_split(mean_error, P)[j]*1.96,
                     fmt=".", ms=20,elinewidth=3,alpha=0.5, linestyle=':')
    ax[0,0].set_ylabel(r"$\langle {P} \rangle$ [mm]")
    ax[0,j].set_title(subplots_text[j])
    ax[0,j].set_ylim(bottom=0)
    ax[0,j].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax[0,j].text(0.03, 0.05, sub_text[j], fontsize=20,transform=ax[0,j].transAxes)

    ax[1,j].errorbar(timesteps, np.array_split(st_dev, P)[j],  color='k',  yerr=np.array_split(st_dev_error, P)[j]*1.96, fmt=".",
                     alpha=0.5,elinewidth=3, ms=20, linestyle=':')
    ax[1,0].set_ylabel(r"$\sigma (P)$ [mm]")
    ax[1,j].set_ylim(bottom=0)
    ax[1,j].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax[1,j].set_xticks(ticks=timesteps, rotation=45, fontsize=10)
    ax[1,j].set_xticklabels(timesteps, rotation=45)
    ax[1,j].text(0.03, 0.05, sub_text[j+P], fontsize=20,transform=ax[1,j].transAxes)
    ax[1,j].set(xlabel='$\Delta t_\mathrm{coal}$ [s]')

plt.subplots_adjust(left=0.05,
                bottom=0.08, 
                right=0.98, 
                top=0.95, 
                wspace=0.22, 
                hspace=0.05)

grid = plt.GridSpec(2, 4)
plt.savefig('PATH TO WHERE TO SAVE FIG.pdf')
