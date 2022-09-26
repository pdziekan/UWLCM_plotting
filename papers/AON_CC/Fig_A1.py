###INIT
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
SD_to_plot_init = [10, 100, 1000, 5000, 10000]
subplots_text = ['tail', 'c multi']

L = int(len(subplots_text))
fig, ax = plt.subplots(2, L, sharex=True)
fig.set_size_inches(19.5, 15.0)

for j in range(len(subplots_text)):        
    ax[0,j].errorbar(SD_to_plot_init, np.array_split(mean_value, L)[j],  color='k',  yerr=np.array_split(mean_error, L)[j]*1.96,
                     fmt=".", ms=20,elinewidth=3,alpha=0.5, linestyle=':')
    ax[0,0].set_ylabel(r"$\langle {P} \rangle$ [mm]")
    ax[0,j].set_xscale('log')
    ax[0,j].set_title(subplots_text[j])
    ax[0,j].text(0.03, 0.05, sub_text[j], fontsize=20,transform=ax[0,j].transAxes)
    ax[0,j].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax[0,j].set_ylim(bottom=0)

    ax[1,j].errorbar(SD_to_plot_init, np.array_split(st_dev, L)[j],  color='k',  yerr=np.array_split(st_dev_error, L)[j]*1.96, fmt=".",
                     alpha=0.5,elinewidth=3, ms=20, linestyle=':')
    ax[1,0].set_ylabel(r"$\sigma (P)$ [mm]")
    ax[1,j].set_xscale('log')
    ax[1,j].set_xlim(1, 10000)
    ax[1,j].set_xlabel('N$\mathrm{_{SD}}$')
    ax[1,j].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax[1,j].get_xaxis().set_minor_formatter(matplotlib.ticker.LogFormatter())
    ax[1,j].set_ylim(bottom=0)
    ax[1,j].text(0.03, 0.05, sub_text[j+L], fontsize=20,transform=ax[1,j].transAxes)


plt.subplots_adjust(left=0.07,
                bottom=0.05, 
                right=0.98, 
                top=0.95, 
                wspace=0.3, 
                hspace=0.05)

grid = plt.GridSpec(2, 4)
plt.savefig('PATH TO WHERE TO SAVE FIG .pdf')
