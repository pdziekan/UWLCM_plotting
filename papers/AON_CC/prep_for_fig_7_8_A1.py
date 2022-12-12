from math import exp, log, sqrt, pi, erf, cos, pow, asin, atan, acos, factorial
import numpy as np
from scipy.stats import moment
from sys import argv
#from matplotlib.ticker import MultipleLocator
import matplotlib as plt
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.colors as mcolors
#from matplotlib.ticker import FormatStrFormatter, LogFormatter, LogLocator, LogFormatterSciNotation, AutoMinorLocator
import glob, os
plt.rcParams.update({'font.size': 20})

############################################################################
####  DO RECZNEGO UZUPELNIENIA#############

label_list = [ '10', '50','100', '1000','10000', '40000', '100000']
SD_list = [ 10, 50, 100, 1000, 10000, 40000, 100000 ] * 4


paths = ['''
           PLEAS PROVIDE PATHS TO RESULTRS FOR EACH SD VALUE WITHIN EACH SCENARIO 

	''']

name = 'Rain_compare'
text_to_legend = ['D', 'LR', 'MR', 'HR']
outfile = 'WRITE A PATHE WHERE TO SAVE DATA/'
width_multiplier = 0.37
##########################################################################
def do_plot(paths_to_files, text_labels, text_to_legend, name):

    if len(text_to_legend) == 1 :
        label = text_labels
    else:
        label = text_labels[0:int(len(paths_to_files)/2)]*len(text_to_legend)
    multi = len(text_to_legend)
    Y = [i+1 for i in range(multi)]
    X = np.repeat(Y, int(len(text_labels)))
    labels = np.repeat(text_to_legend, int(len(text_labels)))
    def read_my_array(file_obj):
        arr_name = file_obj.readline()
        file_obj.readline() # discarded line with size of the array
        line = file_obj.readline()
        line = line.split(" ")
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

    def calc_value(parameter_name, iter_value, paths_to_files):
        dl = len(series_file[iter_value])
        mean_value =[0 for i in range(len(paths_to_files))]
        st_dev = [0 for i in range(len(paths_to_files))]
        mean_error = [0 for i in range(len(paths_to_files))]
        st_dev_error = [0 for i in range(len(paths_to_files))]
        variable = np.zeros((len(series_file[iter_value])))
        for j in range(len(series_file[iter_value])):
            variable[j] = read_my_var(series_file[iter_value][j], str(parameter_name))[-1]
        mean_value[iter_value] = variable.mean(0)
        st_dev[iter_value] = variable.std(0)
        mean_error[iter_value] = variable.std(0)/sqrt(dl)
        st_dev_error[iter_value] = np.power(1/dl * (moment(variable,4) - (dl-3)/(dl-1)*np.power(st_dev[iter_value],4)),1/2)/(2*st_dev[iter_value])
        st_dev_error[iter_value] = np.nan_to_num(st_dev_error[iter_value])
        #Error from this https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
        return mean_value[iter_value], st_dev[iter_value], st_dev_error[iter_value], mean_error[iter_value]

    def sim_time(iter_value):
        for j in range(len(series_file[iter_value])):
            calc_time = read_my_var(series_file[iter_value][j], "position")
        return calc_time

    files = [0 for i in range(len(paths_to_files))]
    series_file = [0 for i in range(len(paths_to_files))]
    for p in range(len(paths_to_files)):
        os.chdir(paths_to_files[p])
        series_file[p] = [open(file_names, "r") for file_names in glob.glob("*.dat")]
        files[p] = glob.glob("*.dat")

    colors = [ 'forestgreen', 'gold', 'forestgreen', 'gold', 'forestgreen', 'gold']
    u_init = np.linspace(-int(len(label)/len(text_to_legend)), int(len(label)/len(text_to_legend)), int(len(label)/len(text_to_legend)))
    u = np.tile(u_init, len(text_to_legend))
    width = 0.20
    colors_list = [ 'thistle', 'orchid','red', 'grey', 'green', 'orange', 'blue', 'yellow', 'forestgreen', 'gold', 'navy', 'crimson' ]
    colors = colors_list[0:len(text_labels)] *len(text_to_legend)
    multi = len(text_labels)

    mean_value, st_dev, st_dev_error, mean_error = [[0 for i in range(len(paths))] for i in range(4)]
    for n in range(len(paths)):
        mean_value[n], st_dev[n], st_dev_error[n], mean_error[n] = calc_value("acc_precip", n, paths_to_files)

    np.save(outfile+'mean_value', mean_value)
    np.save(outfile+'st_dev', st_dev)
    np.save(outfile+'st_dev_error', st_dev_error)
    np.save(outfile+'mean_error', mean_error)
    np.save(outfile+'paths', paths)

do_plot(paths, label_list, text_to_legend, name)
