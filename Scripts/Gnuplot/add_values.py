from __future__ import division

import pylab as pl
import numpy as np
import os,re,sys,getopt
import numpy.ma as MA 

#Arguments
data_file = sys.argv[1] #Folder name: file.gid
var_x     = sys.argv[2] #Variable: variable to add value to
add_val   = sys.argv[3] #Variable: value to add
cols      = sys.argv[4] #Variable: column names


# Open file and read column names and data block 
f  = open(sys.argv[1]) 

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
for (line_count, line) in enumerate(data_block): 
        data[col_names[int(var_x)]][line_count] = float(data[col_names[int(var_x)]][line_count]) + float(add_val)
#return data 

fw = open(sys.argv[1]+'val','w')

for (line_count, line) in enumerate(data_block): 
    items = line.split() 
    for (col_count, col_name) in enumerate(col_names): 
        print>>fw,'  ','{:0=8.07e}'.format(data[col_name][line_count]),
        
    print>>fw,' '

fw.close() 
