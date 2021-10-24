from __future__ import division
import scipy.fftpack

import subprocess
import numpy as np
import os,sys,getopt
import numpy.ma as MA 

#Default Arguments
folder   = ""      #'turek.gid' for example
file_n   = ""      #'turek.sld.tp' for example
cols     = ""      # column names 'Time Vel_x Vel_y Press'
skip     = 0       # Lines to skip in case it has a header 0
var_x    = 0       # variable to plot on Y axis 0 --> Time, default to time
var_y    = ""      # variable to plot on Y axis 1 --> Vel_x, no default
llim     = -1e12   # lower limit to plot on x axis
ulim     =  1e12   # upper limit to plot on x axis
mode     = "multi" # Processing many or single files at once
key      = "yes"   # Print key on plots
doFFT    = "no"    # Do FFT of results
subx     = "no"    # substitute x for linear space from 0 to len(y)
printall = "no"    # print all variables

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:f:x:y:c:s:b:t:m:k:z:p:l:",["directory=","filename=","xaxis=","yaxis=","column_names=","skip=","lowerlimit=","upperlimit=","mode=","key=","fft=","subx=","printall="])
except getopt.GetoptError:
    print 'Wrong options, usage: multignuplot.py -d <directory> -f <file_name> -x <X_axis> -y <Y_axis> -c <column_names> -s <rows_to_skip> -b <lower_limit_to_plot> -t <upper_limit_to_plot> -z <fft> -k <key> -m <mode> -l <subx> -p <printall>'
    sys.exit(2)

for opt, arg in opts:
  if opt == '-h':
    print 'Help: multignuplot.py -d <directory> -f <file_name> -x <X_axis> -y <Y_axis> -c <column_names> -s <rows_to_skip> -b <lower_limit_to_plot> -t <upper_limit_to_plot> -k <key> -z <fft> -m <mode> -l <subx> -p <printall>'
    sys.exit()
  elif opt in ("-d", "--directory"):
    folder     = arg
  elif opt in ("-f", "--inputfile"):
    file_n     = arg
  elif opt in ("-x", "--xaxis"):
    var_x      = arg
  elif opt in ("-y", "--yaxis"):
    var_y      = arg
  elif opt in ("-c", "--colnames"):
    cols       = arg
  elif opt in ("-s", "--skip"):
    skip       = arg
  elif opt in ("-b", "--lowerlimit"):
    llim      = float(arg)
  elif opt in ("-t", "--upperlimit"):
    ulim      = float(arg)
  elif opt in ("-m", "--mode"):
    mode      = arg
  elif opt in ("-k", "--key"):
    key       = arg
  elif opt in ("-z", "--fft"):
    doFFT     = arg
  elif opt in ("-l", "--subx"):
    subx      = arg
  elif opt in ("-p", "--printall"):
    printall  = arg

#Array of result files
f_names = [] 
col_names = [] 
#Array of names to label results
k_names = [] #not the same as above as some files could be empty
# Create a data dictionary, containing 
# a list of values for each variable 
data    = {} 
#Array of root names to label results
r_names = {} #Root folder names used to change name if necessary
change  = [] #If a result file is in subfolder we adopt the root name

def processFolders(mode,folder,cols,f_names,col_names):

    global change

    if(mode == "multi"):
        os.chdir(folder)
    
        col_names.extend(cols.split())
        data = {} 
    
    
        #Open file and read column names and data block 
        #We explore all this folder in search for the file
        for root, dirs, files in os.walk(".", topdown=False):
            for name in dirs:
    
                f = os.path.join(root, name + "/" + file_n)
                ig = os.path.join(root, name + "/.ignore")
                ignore = os.path.isfile(ig)
    
                #We found the file, then we open it
                if os.path.exists(f) and not ignore:
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
                print 'develop name change better!'
                data[r_names[str(change)]]=data.pop(change)

        #print data
        return data
    
    elif(mode == "single"):
        # Open file and read column names and data block 
        os.chdir(folder)
        f = open(file_n) 
    
        #col_names = f.readline().split() 
        for i in range(int(skip)):  
            f.readline()  
        data_block = f.readlines() 
        f.close() 
    
        col_names.extend(cols.split())
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
                value = items[col_count] #Time is always col 0 
                data[col_name][line_count] = value 
    
        i=0
        if(printall == "no"):
            f_names.append(col_names[int(var_y)])
        else:
            for i  in range(0,len(col_names)): 
                f_names.append(col_names[int(i)])

        return data
    
    else:
        print 'Wrong processing mode: multi or single available'
        sys.exit(2)

def createPrintFolder(plotdir):
    if os.path.exists(plotdir):
        os.chdir(plotdir)
    else:
        os.mkdir(plotdir)
        os.chdir(plotdir)

def setGnuplotParams(var_y,col_name,proc):
    #Plotting using gnuplot
    proc.stdin.write("set macros\n")
    #proc.stdin.write("load 'dark2.pal' \n")
    #proc.stdin.write("set terminal wxt size 1024,768 font 'Verdana,10' persist \n")
    #proc.stdin.write("set terminal postscript eps enhanced size 12cm,9.7cm color blacktext 'Helvetica, 18' \n")
    proc.stdin.write("set terminal pdfcairo size 12cm,9.7cm \n")

    nam_aux = col_name
    if doFFT == "yes":
        nam_aux = nam_aux + '_fft'
        proc.stdin.write("set logscale x\n") 
        proc.stdin.write("set format x '10^{%01T}'\n") 
        proc.stdin.write("set logscale y\n") 
        proc.stdin.write("set format y '10^{%01T}'\n") 
    if  subx == "yes":
        proc.stdin.write("set logscale y\n") 
        proc.stdin.write("set format y '10^{%01T}'\n") 
    
    #proc.stdin.write("set output 'plot_%s.eps' \n" % nam_aux)
    proc.stdin.write("set output 'plot_%s.pdf' \n" % nam_aux)
    #proc.stdin.write("set terminal pngcairo size 1024,768 enhanced font 'Helvetica, 12' \n")
    #proc.stdin.write("set terminal wxt size 350,262 enhanced font 'Verdana, 12' persist \n")
    
    proc.stdin.write("unset key \n")
    proc.stdin.write("set tics\n")
    #proc.stdin.write("set xrange [0:10] \n") 
    #proc.stdin.write("set yrange [-0.1:0.1] \n")
    
    #proc.stdin.write("set lmargin at screen 0.0; set rmargin at screen 1.9 \n")
    #proc.stdin.write("set tmargin at screen 0.00; set bmargin at screen 1.0 \n")
    proc.stdin.write("set bmargin at screen 0.22 \n")

def setXlabel(var_y,proc):

    lab_x = "Time(s)"
    if doFFT == "yes":
        lab_x = "Frequency (Hz)"
    if subx == "yes":
        lab_x = col_names[int(var_x)]
    proc.stdin.write("set xlabel '%s' \n"%lab_x )
    proc.stdin.write("set ylabel '%s' \n" % col_names[int(var_y)])
    proc.stdin.write("set multiplot layout 1,1 \n")

#Label for graph
#proc.stdin.write("set label 1 '{/Symbol h}_1' at graph 1.0,1.0 \n")

def setAxisData(color,count,mode,data,name,col_x,col_y):

    global change

    command = ""
    #print 'f_names = %s'%(f_names)
    #for name  in f_names: 
    
    if change:
        if name == change:
            name = r_names[str(change)]
    
    if(mode == "multi"):
        x   = data[str(name)][col_x]
        y   = data[str(name)][col_y]
        f = 'plot.'+ (name.replace("./","")).split("/", 1)[0]
    elif(mode == "single"):
        x   = data[col_x]
        y   = data[col_y]
        f = 'plot.'+ name[count]
    
    x_=[]
    y_=[]
    if (subx=="yes"):
        x     = np.linspace(0, len(y),len(y))
    
    for j in range(1,len(x)):
        if (subx=="yes"):
            if x[j] >=llim and x[j] <= ulim:
                x_.append(j)
                y_.append(y[j])
        else:
            if x[j] >=llim and x[j] <= ulim:
                x_.append(x[j])
                y_.append(y[j])
    
    if len(x) > 1:
        if doFFT == "yes":
            sample  = len(x_)
            spacing = x_[len(x_)-2]-x_[len(x_)-3]
            x_     = np.linspace(0.0, 1.0/(2.0*spacing), sample/2)
            y_aux  = abs(scipy.fftpack.fft(y_))
            y_     = 2.0/sample * np.abs(y_aux[:sample//2])
            f = f+'_fft'
    
    temp = open(f,"w+")
    for j in range(1,len(x_)):
        temp.write("%s %s \n" % (x_[j],y_[j]))
    temp.close() 
    
    #proc.stdin.write("plot '%s' using 1:2  with l ls '%s' lw 1; " % (f,i))
    if os.stat(f).st_size != 0:
        k_names.append(name)

        if(mode == "single"):
            aux     = "'%s' using 1:2 with l ls %i lw 1" % (f,color+1)
        elif(mode == "multi"):
            aux     = "'%s' using 1:2 with lp ls %i lw 1 pi 50 pt %i ps 1.0" % (f,color+1,2*(color+1))

        command = ",".join((aux,command))

    return command

def plotData(command,proc):
    #Take out the last comma
    command = command[:-1]
    
    proc.stdin.write("plot %s " %command )
    proc.stdin.write(" \n" )

def PrintPlotLabels(key,proc,count):
    global change
    if (key=="yes"):
        ### Key plot
        proc.stdin.write("set key horizontal bottom center\n")
        proc.stdin.write("set border 0\n")
        proc.stdin.write("unset label\n")
        proc.stdin.write("unset tics\n")
        proc.stdin.write("unset xlabel\n")
        proc.stdin.write("unset ylabel\n")
        proc.stdin.write("set yrange [0.1:0.2]\n")
        proc.stdin.write("set xrange [0.1:0.2]\n")
        proc.stdin.write("set lmargin at screen 0.02; set rmargin at screen 1.95 \n")
        proc.stdin.write("set tmargin at screen 0.40; set bmargin at screen 0.035 \n")
        
        i = 1
        command = ""
        for name  in k_names: 
    
            if change:
                if name == change:
                    name = r_names[str(change)]
    
            if isinstance(name, basestring):
                newname = (name.replace("./","")).split("/", 1)[0]
            else:
                newname = (name[count-1].replace("./","")).split("/", 1)[0]


            #aux     = "1 t '%s' noenhanced with lp ls %i" % (newname,i)
            aux     = "1 t '%s' noenhanced with lp ls %i pt %i" % (newname,i,2*i)
            command = ",".join((aux,command))
            i = i + 1
            
        command = command[:-1]
        proc.stdin.write("plot %s " %command )

def plot():

    global k_names

    if(printall == "yes"):

        if(mode == "multi"):

            for j  in range(1,len(col_names)): 
                proc = subprocess.Popen(['gnuplot','-p'],shell=True, stdin=subprocess.PIPE,preexec_fn=os.setsid)
                setGnuplotParams(j,col_names[j],proc)
                setXlabel(j,proc)
                k_names =[]

                command = ""
                for i  in range(0,len(f_names)): 
                    #print 'name :',f_names[i],'; y col :',col_names[j]
                    aux     = setAxisData(i,i,mode,data,f_names[i],col_names[int(var_x)],col_names[j])
                    command = "".join((aux,command))

                plotData(command,proc)
                PrintPlotLabels(key,proc,i)
                #os.killpg(os.getpgid(proc.pid), signal.SIGTERM)  # Send the signal to all the process groups
                proc.kill()

        else:

            for i  in range(1,len(col_names)): 
                proc = subprocess.Popen(['gnuplot','-p'],shell=True, stdin=subprocess.PIPE,preexec_fn=os.setsid)
                setGnuplotParams(i,col_names[i],proc)
                setXlabel(i,proc)
                k_names =[]
                command = ""
                command = setAxisData(0,i,mode,data,f_names,col_names[int(var_x)],col_names[i])
                plotData(command,proc)
                PrintPlotLabels(key,proc,i)
                #os.killpg(os.getpgid(proc.pid), signal.SIGTERM)  # Send the signal to all the process groups
                proc.kill()
    else:

        proc = subprocess.Popen(['gnuplot','-p'],shell=True, stdin=subprocess.PIPE)
        setGnuplotParams(var_y,col_names[int(var_y)],proc)
        setXlabel(var_y,proc)
        command = ""
        aux= ""

        if(mode == "multi"):
            for i  in range(0,len(f_names)): 
                aux= setAxisData(i,i,mode,data,f_names[i],col_names[int(var_x)],col_names[int(var_y)])
                command = "".join((aux,command))
        else:
            command = setAxisData(0,0,mode,data,f_names,col_names[int(var_x)],col_names[int(var_y)])

        plotData(command,proc)
        PrintPlotLabels(key,proc,0)
        proc.kill()

#----------------------------------Main------------------------------------ 
data=processFolders(mode,folder,cols,f_names,col_names)
plotdir = "gnuplots"
createPrintFolder(plotdir)
plot()

