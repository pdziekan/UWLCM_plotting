from math import exp, log, sqrt, pi, erf, cos, pow, asin, atan, acos, factorial
import numpy as np
from sys import argv
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter, LogFormatter, LogLocator, LogFormatterSciNotation, AutoMinorLocator
import glob, os
plt.rcParams.update({'font.size': 18})

def multiplot(name, paths, label, outfile, sd):
    def read_my_array(file_obj):
        arr_name = file_obj.readline()
        file_obj.readline() # discarded line with size of the array
        line = file_obj.readline()
        line = line.split(" ")
        print(len(line), line)
        del line[0]
        del line[len(line)-1]
        arr = list(map(float,line))
        return np.array(arr), arr_name

    def read_my_var(file_obj, var_name):
        file_obj.seek(0)
        while True:
            arr, name = read_my_array(file_obj)
            if(str(name).strip() == str(var_name).strip()):
                break
        return arr

    def calc_value(parameter_name, iter_value, paths):
        dl = len(parameter_name)
        mean_value =[0 for i in range(len(paths))]
        st_dev = [0 for i in range(len(paths))]
        variable = np.zeros((len(series_file[iter_value]),len(read_my_var(series_file[iter_value][0], str(parameter_name)))))
        for j in range(len(series_file[iter_value])):
            variable[j] = read_my_var(series_file[iter_value][j], str(parameter_name))
        mean_value[iter_value] = variable.mean(0)
        st_dev[iter_value] = variable.std(0)
        return mean_value[iter_value], st_dev[iter_value]

    def sim_time(iter_value):
        for j in range(len(series_file[iter_value])):
            calc_time = read_my_var(series_file[iter_value][j], "position")
        return calc_time

    files = [0 for i in range(len(paths))]
    series_file = [0 for i in range(len(paths))]
    for p in range(len(paths)):
        os.chdir(paths[p])
        series_file[p] = [open(file_names, "r") for file_names in glob.glob("*.dat")]
        files[p] = glob.glob("*.dat")

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    colors = ['red', 'green','blue', 'orange', 'green', 'orange', 'blue', 'yellow', 'forestgreen', 'gold', 'navy', 'crimson']

    fig.set_size_inches(14.5, 20.5)
    axis1 = plt.subplot(511)
    for p in range(len(paths)):
        cth = calc_value("cloud_top_height", p, paths)[0]
        time = sim_time(p)
        plt.plot(time, cth,  color=colors[p], linewidth=3, label=(label[p]))
        plt.ylim((0, 5000))
        plt.ylabel("CTH [m]")
        plt.text(100,5000*0.9,"(a)",fontsize = 22)
        plt.legend(title=f"Average of different scenarios realizations with {sd} SD per gridbox :", loc='upper center',
                bbox_to_anchor=(0.5, 1.55), ncol=4,frameon=0 )
        plt.tick_params('x', labelbottom=False)
    axis2 = plt.subplot(512, sharex=axis1)
    for p in range(len(paths)):
        cc = calc_value("cloud_cover_rico", p, paths)[0]
        time = sim_time(p)
        plt.plot(time, cc,  color=colors[p], linewidth=3, label=label[p])
        plt.ylabel("Cloud Cover [-]")
        plt.ylim((0, 0.3))
        plt.text(100,0.3*0.9,"(b)",fontsize = 22)
        plt.tick_params('x', labelbottom=False)
    axis3 = plt.subplot(513, sharex=axis1)
    for p in range(len(paths)):
        cwp = calc_value("cwp", p, paths)[0]
        time = sim_time(p)
        plt.plot(time, cwp,  color=colors[p], linewidth=3, label=(label[p]))
        plt.ylabel("CWP [g/m^2]")
        plt.ylim((0, 200))
        plt.text(100,200*0.9,"(c)",fontsize = 22)
        plt.tick_params('x', labelbottom=False)
    axis4 = plt.subplot(514, sharex=axis1)
    for p in range(len(paths)):
        rwp = calc_value("rwp", p, paths)[0]
        time = sim_time(p)
        plt.plot(time, rwp,  color=colors[p], linewidth=3, label=(label[p]))
        plt.ylabel("RWP [g/m^2]")
        plt.xlim((0, 10850))
        plt.ylim((0, 220))
        plt.tick_params('x', labelbottom=False)
        plt.text(100,220*0.9,"(d)",fontsize = 22)
    axis5 = plt.subplot(515)
    for p in range(len(paths)):
        sur_p = calc_value("surf_precip", p, paths)[0]
        cc = calc_value("cloud_cover_rico", p, paths)[0]
        time = sim_time(p)
        plt.plot(time, sur_p/cc/24,  color=colors[p], linewidth=4, label=(label[p]))
        plt.xlabel("time [s]")
        plt.ylabel("Precipitation [mm/h] ")
        plt.xlim((0, 10850))
        plt.yscale('log')
        plt.ylim((10**-2, 2*10))
        plt.text(100,2*10*0.5,"(e)",fontsize = 22)
    fig.tight_layout()
    fig.savefig(outfile + 'CC_properties_'+name+'.pdf')


SD = [10, 50, 100, 1000, 10000, 40000, 100000]

for sd in SD:
    name_to_plot = f'Fig_6_SD{sd}'
    paths_to_plot = [f'REPLACE WITH THE PATH TO D SCENARIO RESULTS FOR {sd}',
            f'REPLACE WITH THE PATH TO LR SCENARIO RESULTS FOR {sd} ',
            f'REPLACE WITH THE PATH TO MR SCENARIO RESULTS FOR {sd}',
            f'REPLACE WITH THE PATH TO HR SCENARIO RESULTS FOR {sd}']
    label_to_plot = ['D','LR','MR','HR']
    outfile_to_plot = 'REPLACE WITH THE PATH WHERE TO SAVE PLOT'

    multiplot(name_to_plot, paths_to_plot, label_to_plot, outfile_to_plot, sd)
