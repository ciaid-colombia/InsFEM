from __future__ import division

import pylab as pl
import numpy as np
import os,re,sys,getopt
import numpy.ma as MA 

#Arguments
data_file = sys.argv[1] #Folder name: file.gid
cols      = sys.argv[2] #Variable: column names
var_x     = sys.argv[3] #Variable: variable to plot on Y axis
var_y     = sys.argv[4] #Variable: variable to plot on Y axis

#var_x     = ''
#var_y     = ''
#cols      = ''
#if len(sys.argv) == 1:
# print 'plot_tp.py -f <file_name> -x <X_axis> -y <Y_axis> -c <column_names>'
# sys.exit()
#try:
#  opts, args = getopt.getopt(sys.argv[1:],"hf:x:y:c:",["filename=","xaxis=","yaxis=","column_names"])
#except getopt.GetoptError:
#  print 'plot_tp.py -f <file_name> -x <X_axis> -y <Y_axis> -c <column_names>'
#  sys.exit(2)
#for opt, arg in opts:
#  if opt == '-h':
#    print 'plot_tp.py -f <file_name> -x <X_axis> -y <Y_axis> -c <column_names>'
#    sys.exit()
#  elif opt in ("-f", "--inputfile"):
#    f     = arg
#  elif opt in ("-x", "--xaxis"):
#    var_x = arg
#  elif opt in ("-y", "--yaxis"):
#    var_y = arg
#  elif opt in ("-c", "--colnames"):
#    cols  = arg
#
#print(len(sys.argv))
#print(var_x)
#print(var_y)
#print(cols)

# Open file and read column names and data block 
f = open(sys.argv[1]) 

#col_names = f.readline().split() 
data_block = f.readlines() 
f.close() 

col_names = cols.split()
# Create a data dictionary, containing 
# a list of values for each variable 
data = {} 

# Add an entry to the dictionary for each column 
for col_name in col_names: 
    data[col_name] = MA.zeros(len(data_block), 'f',  
            fill_value = -999.999) 

# Loop through each value: append to each column 
for (line_count, line) in enumerate(data_block): 
    items = line.split() 
    for (col_count, col_name) in enumerate(col_names): 
        value = items[col_count] 
        data[col_name][line_count] = value 
#return data 

dt   = data[col_names[int(var_x)]]
Y    = data[col_names[int(var_y)]]

#Plotting
from pylab import plot
fig = pl.figure()
ax = fig.add_subplot(111)
ax.ticklabel_format(axis='y', style='sci')
ax.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))
#ax.set_yscale('log')
pl.plot(dt,Y,color="red", linewidth=1.0, linestyle="-",label = col_names[int(var_y)])
#pl.xlim(0,120)
pl.xlabel('Time (s)')
pl.ylabel(col_names[int(var_y)])
pl.grid()
#Two options, show print on screen, savefig saves the file to disk
pl.show()
#pl.savefig(suffix + 'fft.png',dpi=150)

