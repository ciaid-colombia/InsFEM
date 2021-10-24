from __future__ import division

import pylab as pl
import numpy as np
import os,re,sys,getopt
import numpy.ma as MA 
import random as rm
from pylab import plot

#Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#RGB color; the keyword argument name must be a standard mpl colormap name
def get_cmap(n, name='gnuplot'):
    return pl.cm.get_cmap(name, n)


#Arguments
folder = sys.argv[1] #'turek.sld.tp' for example
file_n = sys.argv[2] #'turek.sld.tp' for example
cols   = sys.argv[3] # column names 'Time Vel_x Vel_y Press'
skip   = sys.argv[4] # Lines to skip in case it has a header 0
var_x  = sys.argv[5] # variable to plot on Y axis 0 --> Time
var_y  = sys.argv[6] # variable to plot on Y axis 1 --> Vel_x

os.chdir(folder)

col_names = cols.split()

# Create a data dictionary, containing 
# a list of values for each variable 
data = {} 
#Array of names to label results
f_names = [] 
#Array of root names to label results
r_names = {}
change = []

# Open file and read column names and data block 
#We explore all this folder in search for the file
for root, dirs, files in os.walk(".", topdown=False):
    for name in dirs:

        f = os.path.join(root, name + "/" + file_n)

        #We found the file, then we open it
        if os.path.exists(f):
            f_tp = open(f,"r") 

            #We check to see if name is already mapped
            if data.has_key(name):
                #We replace name for root name and change the original
                change=str(name)
                name = root

            r_names.update({name:root})
            f_names.append(name)

            #Skip header lines
            for i in range(int(skip)):  
                f_tp.readline()  

            data_block = f_tp.readlines() 
            f_tp.close() 

            # Add an entry to the dictionary for each column of each case
            aux_data={name:{}}
            data.update(aux_data)
            data[str(name)].fromkeys(col_names)
            for col_name in col_names: 
                data[str(name)][col_name] = MA.zeros(len(data_block), 'f',fill_value = -999.999)

                # Loop through each value: append to each column 
            for (line_count, line) in enumerate(data_block): 
                items = line.split() 

                for (col_count, col_name) in enumerate(col_names): 
                    value = items[col_count] 
                    data[str(name)][col_name][line_count]=value 

if change:
    data[r_names[str(change)]]=data.pop(change)

i=1
#Plotting
fig = pl.figure()
ax = fig.add_subplot(111)
ax.ticklabel_format(axis='y', style='sci')
#Get random color array
cmap = get_cmap(len(f_names)*10)
#ax.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))
#ax.set_yscale('log')
pl.hold(True)
pl.xlabel('Time (s)')
pl.grid()

for name  in f_names: 

    if change:
        if name == change:
            name = r_names[str(change)]

    x   = data[str(name)][col_names[int(var_x)]]
    y   = data[str(name)][col_names[int(var_y)]]

    #pl.plot(x,y,color=cmap(i*rm.randint(1,10)), linewidth=1.0, linestyle="-",label = str(name))
    pl.plot(x,y,color=cmap(i*8), linewidth=1.0, linestyle="-",label = str(name))
    pl.legend(loc='best').draggable()
    i=i+1
    
    
#Two options, show print on screen, savefig saves the file to disk
pl.show()
#pl.savefig(file_n+ 'fft.png',dpi=150)


