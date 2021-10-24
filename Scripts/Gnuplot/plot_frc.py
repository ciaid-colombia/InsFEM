from __future__ import division

import pylab as pl
import numpy as np
import os,re,sys
import numpy.ma as MA 

#Arguments
#data_file = sys.argv[1] #Folder name: file.gid
var       = sys.argv[2] #Variable: variable to plot on Y axis
suffix    = sys.argv[3] #Variable: name of variable to label Y axis

# Open file and read column names and data block 
f = open(sys.argv[1]) 

# Ignore header 
for i in range(1):  
    f.readline()  
col_names = f.readline().split() 
data_block = f.readlines() 
f.close() 

col_names = col_names[1:5]
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

dt   = data["Time"]
Y    = data[var]

#Plotting
from pylab import plot
fig = pl.figure()
ax = fig.add_subplot(111)
ax.ticklabel_format(axis='y', style='sci')
ax.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))
#ax.set_yscale('log')
pl.plot(dt,Y,color="red", linewidth=1.0, linestyle="-",label = suffix)
#pl.xlim(0,120)
pl.xlabel('Time (s)')
pl.ylabel(suffix)
pl.grid()
#Two options, show print on screen, savefig saves the file to disk
pl.show()
#pl.savefig(suffix + 'fft.png',dpi=150)

