#!/usr/bin/env python

import datetime
import math
import numpy as np
try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk
try:
    import tkFileDialog as fd
except ImportError:
    import tk.filedialog as fd
try:
    import tkMessageBox as mb
except ImportError:
    import messagebox as mb
import sys
from StringIO import StringIO
import struct
import os
import csv
import time as time_module
#import f2py2e

##Change the current directory to the one that the script is located
##REFERENCE: Stack Overflow, Eli Courtwright

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

##Collect the date/time to label the subdirectory

date = datetime.datetime.now().strftime("%m_%d_%Y_%H_%M_%S")

os.makedirs(date)
os.chdir(date)
os.system('cp ../new_align.f ./')
os.system('cp ../Monte_Carlo.f ./')


##Compile the fortran scripts.  This is done early in the script in order to prevent the
##user from running most of the script before a potential failure.

os.system("gfortran -C new_align.f -o new_align.out")

##The above command can fail, if the user does not have command line developer tools installed.
##If this happens, it will prompt the user to install it, then quit the program.

if os.path.isfile("new_align.out"):
	pass
else:
	sys.exit("The script was unable to compile the fortran code.  If your computer prompts you to install Command Line Developer Tools, please do so and then restart the program.  If you have not already installed gfortran, please do so as well.")

os.system("gfortran -llapack -lblas Monte_Carlo.f -o Monte_Carlo.out")

#np.lookfor('f2py')

#np.f2py -c -m new_align new_align.f

#np.f2py.run_main('new_align.f')
#with open("new_align.f90") as sourcefile:
#    sourcecode = sourcefile.read()
#np.f2py.compile(sourcecode, modulename='new_align',extra_args = '--fcompiler=gfortran')
#with open("Monte_Carlo.f") as sourcefile:
#	sourcecode = sourcefile.read()
#np.f2py.compile(sourcecode, modulename='Monte_Carlo',extra_args = '--fcompiler=gfortran')

#import new_align
#import Monte_Carlo
#f2py -m -c new_align 'new_align.f'
#f2py -c -lblas -llapack Monte_Carlo.out Monte_Carlo.f


##Define fcsextract function



########################## fcsextract #############################

def fcsextract(filename):
    """
    Attempts to parse an FCS (flow cytometry standard) file

    Parameters: filename
        filename: path to the FCS file

    Returns: (vars,events)
    	vars: a dictionary with the KEY/VALUE pairs found in the HEADER
    	this includes the standard '$ABC' style FCS variable as well as any 
    	custom variables added to the header by the machine or operator
	
    	events: an [N x D] matrix of the data (as a Python list of lists)
    	i.e. events[99][2] would be the value at the 3rd dimension
    	of the 100th event
    """
    fcs_file_name = filename

    fcs = open(fcs_file_name,'rb')
    header = fcs.read(58)
    version = header[0:6].strip()
    text_start = int(header[10:18].strip())
    text_end = int(header[18:26].strip())
    data_start = int(header[26:34].strip())
    data_end = int(header[34:42].strip())
    analysis_start = int(header[42:50].strip())
    analysis_end = int(header[50:58].strip())

    print "Parsing TEXT segment"
    # read TEXT portion
    fcs.seek(text_start)
    delimeter = fcs.read(1)
    # First byte of the text portion defines the delimeter
    print "delimeter:",delimeter
    text = fcs.read(text_end-text_start+1)

    #Variables in TEXT poriton are stored "key/value/key/value/key/value"
    keyvalarray = text.split(delimeter)
    fcs_vars = {}
    fcs_var_list = []
    # Iterate over every 2 consecutive elements of the array
    for k,v in zip(keyvalarray[::2],keyvalarray[1::2]):
        fcs_vars[k] = v
        fcs_var_list.append((k,v)) # Keep a list around so we can print them in order

    #from pprint import pprint; pprint(fcs_var_list)
    if data_start == 0 and data_end == 0:
        data_start = int(fcs_vars['$DATASTART'])
        data_end = int(fcs_vars['$DATAEND'])

    num_dims = int(fcs_vars['$PAR'])
    print "Number of dimensions:",num_dims

    num_events = int(fcs_vars['$TOT'])
    print "Number of events:",num_events

    # Read DATA portion
    fcs.seek(data_start)
    #print "# of Data bytes",data_end-data_start+1
    data = fcs.read(data_end-data_start+1)

    # Determine data format
    datatype = fcs_vars['$DATATYPE']
    if datatype == 'F':
        datatype = 'f' # set proper data mode for struct module
        print "Data stored as single-precision (32-bit) floating point numbers"
    elif datatype == 'D':
        datatype = 'd' # set proper data mode for struct module
        print "Data stored as double-precision (64-bit) floating point numbers"
    else:
        assert False,"Error: Unrecognized $DATATYPE '%s'" % datatype
    
    # Determine endianess
    endian = fcs_vars['$BYTEORD']
    if endian == "4,3,2,1":
        endian = ">" # set proper data mode for struct module
        print "Big endian data format"
    elif endian == "1,2,3,4":
        print "Little endian data format"
        endian = "<" # set proper data mode for struct module
    else:
        assert False,"Error: This script can only read data encoded with $BYTEORD = 1,2,3,4 or 4,3,2,1"

    # Put data in StringIO so we can read bytes like a file    
    data = StringIO(data)

    print "Parsing DATA segment"
    # Create format string based on endianeness and the specified data type
    format = endian + str(num_dims) + datatype
    datasize = struct.calcsize(format)
    print "Data format:",format
    print "Data size:",datasize
    events = []
    # Read and unpack all the events from the data
    for e in range(num_events):
        event = struct.unpack(format,data.read(datasize))
        events.append(event)
    
    fcs.close()
    return fcs_vars,events
    
def writefcs(fcs_vars,events,fcs_file_name,delimiter=","):
    """
    Outputs FCS variables and data to files
    
    fcs_vars: the dictionary of key/value pairs from HEADER
    events: [N x D] matrix (list of lists) of event data in row-major form
    fcs_file_name: prefix for the output files
    delimiter: specifies separator between values in ASCII file output
        Generates a binary file if None
        
    Creates 3 files
    a) HEADER: fcs_file_name.txt
        the HEADER key/value pairs
    b) DATA: fcs_file_name.csv (or .bin for binary file)
        the raw data, one event per line
    c) INFO: fcs_file_name.info
        list of the dimension names and long-names ($PkN and $PkS)
    """
    num_dims = len(events[0])
    num_events = len(events)

    if delimiter is None:
        # Creates a binary file
        # First 4 bytes are an integer with the number of events
        # Next 4 bytes are an integer with the number of dimensions
        # Rest of the file is consecutive 32-bit floating point numbers
        # Data is stored such that consecutive floats are from the same event 
        # (i.e. an N x D matrix in row-major format)
        bin_file_name = fcs_file_name[:-4]+".bin"
        bin_file = open(bin_file_name,"wb")
        print "Writing DATA output file:",bin_file_name
        bin_file.write(struct.pack("i",num_events))
        bin_file.write(struct.pack("i",num_dims))
        format = "%df" % num_dims
        for row in events:
            data = [float(x) for x in row]
            bin_file.write(struct.pack(format,*data))
    else:
        csv_file_name = fcs_file_name[:-4]+".csv"
        csv_file = open(csv_file_name,'w')
        print "Writing DATA output file:",csv_file_name
        format = delimiter.join(["%0.2f"]*num_dims)
        for event in events:
            csv_file.write(format % event + "\n")
        csv_file.close()

    txt_file_name = fcs_file_name[:-4]+".txt"
    txt_file = open(txt_file_name,'w')
    print "Writing TEXT output file:",txt_file_name
    for k,v in fcs_vars.items():
        txt_file.write("%s,%s\n" % (k,v))
    txt_file.close()

    info_file_name = fcs_file_name[:-4]+".info"
    print "Writing INFO output file:",info_file_name
    info_file = open(info_file_name,'w')
    for i in range(num_dims):
        info_file.write("%s\t%s\n" % (fcs_vars["$P%dN"%(i+1)],fcs_vars.get("$P%dS"%(i+1),fcs_vars["$P%dN"%(i+1)])))     

########################### End fcsextract ####################################

##Execute function for dialogs.  Includes data time inputs, file inputs, and protein selection

############### THE SCRIPT FOR python_dialogs() IS PASTED BELOW ###########################

t = []
f = []

## Prompt the user for times and corresponding files                                                     
class TimeDialog(tk.Frame):

    def __init__(self,master=None):
        tk.Frame.__init__(self,master)
#        mb.showinfo("Time Selection","Input two times for which data was recorded on the following screen.")
        self.get_time()
        

    def get_time(self):
        self.title = tk.Label(text="Input two times for which data was recorded.")
        self.prompt1 = tk.Label(text="First time point:")
        self.prompt2 = tk.Label(text="Second time point:")
#        w.pack(side=tk.LEFT)
        self.entry1 = tk.Entry()
        self.entry2 = tk.Entry()
 
        self.button = tk.Button(text="OK",command=self.ok)
        self.exit_button = tk.Button(text="Quit",command=self.quit_program)

        self.title.grid(row=0,column=0)
        self.prompt1.grid(row=1,column=0)
        self.entry1.grid(row=1,column=1)
        self.prompt2.grid(row=2,column=0)
        self.entry2.grid(row=2,column=1)
        self.button.grid(row=3,column=1)
        self.exit_button.grid(row=3,column=0)

    def ok(self):
        t.append(int(self.entry1.get()))
        t.append(int(self.entry2.get()))

        root.destroy()
        
    def quit_program(self):
    	root.destroy()
    	sys.exit()

root = tk.Tk()
app = TimeDialog()
app.mainloop()
                                                                   
i=0
fid1=open("time.txt","w")
fid2=open("filename.txt","w")

for time in t:
    root = tk.Tk()
    title_str = "Select file for time t="+str(time)+" on the following screen"
    mb.showinfo("File Selection",title_str)
    f.append(fd.askopenfilename(title=title_str))
    if not f[-1]:
    	sys.exit()
    root.destroy()
    fid1.write(str(time))
    fid1.write("\n")
    fid2.write(f[i])
    fid2.write("\n")
    i=i+1

fid1.write("\b")

fid1.flush()
fid1.close
fid2.close

for filename in f:
    (some_vars,some_events) = fcsextract(filename)
    writefcs(some_vars,some_events,filename,delimiter=" ")
    csv_file = filename.replace(".fcs",".csv")
    txt_file = filename.replace(".fcs",".txt")
    os.rename(csv_file,txt_file)

info_file = f[0].replace(".fcs",".info")
#fid=open("info_file_name.txt","w")
#fid.write(info_file)
#fid.close

fid = open(info_file,"r")
lines = fid.readlines()
proteins = []
for line in lines:
    proteins.append(line.split()[1])
fid.close()
del proteins[0]
del proteins[-1]
del proteins[0]
del proteins[-1]

# Sort the proteins in alphabetical order
proteins2 = sorted(proteins, key=str.lower)

status = []
checked = []
final_list = []

class ListDialog(tk.Frame):

    def __init__(self,root):
        tk.Frame.__init__(self,root)
#        mb.showinfo("Protein Selection","Select the proteins to be analyzed on the following screen.")
        self.root = root
        self.scrollbar = tk.Scrollbar(self,orient=tk.VERTICAL)
        self.text = tk.Text(self,width=100,height=400,yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.text.yview)
        self.scrollbar.pack(side=tk.RIGHT,fill="y")
        self.text.pack(side=tk.LEFT,fill="both",expand=True)
        self.exit_button = tk.Button(root,text="Quit Program",command=self.exit_program)
        self.exit_button.pack()
        self.b = tk.Button(root,text="Done",command=self.quit)
        self.b.pack()


        i=0
        for protein in proteins2:
            status.append(tk.IntVar(value=0))
            checkbox = tk.Checkbutton(self,text=protein,variable=status[i])
            
            self.text.window_create("end",window=checkbox)
            self.text.insert("end","\n")
            i=i+1

    def quit(self):
        for i in status:
            checked.append(i.get())
        for i,boo in enumerate(checked):
            if boo==1:
                final_list.append(proteins2[i])
                
        fid=open("protein_list.txt","w")
        for protein in final_list:
            fid.write("%s " % protein)

        root.destroy()
        
    def exit_program(self):
    	root.destroy()
    	sys.exit()

root = tk.Tk()
ListDialog(root).pack(side="top",fill="both",expand=True)
root.mainloop()

proteins = final_list

############################# END python_dialogs() ##################################



##Declare some constants.  This can be updated later to incorporate more robust options.

file_format = '.fcs'
num_steps = 1
tmpoint = 2


##Draw the network between proteins describing their interactions



####################### Script for draw_network.py #################################

##Initialize some variables.  The height and width can be altered if necessary
root = tk.Tk()
canvas_width = root.winfo_screenwidth()/1.5
canvas_height = root.winfo_screenheight()/1.5
root.destroy()
#canvas_width = 800
#canvas_height = 800
center = [canvas_width/2,canvas_height/2]

collection = proteins
n = len(collection)

connection_data = []

##Much of this code is borrowed from the user "Marcin" on Stack Overflow, on their response to the question "drawing rectangle using mouse events in Tkinter"

class NetworkApp(tk.Frame):
    def __init__(self,root2):
        tk.Frame.__init__(self,root2)
        self.root2 = root2
        mb.showinfo("Network Drawing","Draw the connections between selected proteins to be analyzed on the following screen.")
        self.x = self.y = 0
        self.canvas = tk.Canvas(self, width=canvas_width, height=canvas_height)
        self.canvas.pack(side="top", fill="both", expand=True)
        self.canvas.tag_bind('protein',"<ButtonPress-1>", self.on_button_press)
        self.canvas.bind("<B1-Motion>", self.on_move_press)
        self.canvas.tag_bind('protein',"<ButtonRelease-1>", self.on_button_release)
##Create a button to exit the network creation screen
        quit_button = tk.Button(self, text = 'Finish', command=self.quit)
        quit_button.configure(width=10)
        button_window = self.canvas.create_window((canvas_width-150),(canvas_height-50),window=quit_button)
##Create a button to remove the last line drawn
        self.delete_button = tk.Button(self, text = 'Delete last line', command=self.remove_last)
        self.delete_button.configure(width=15,state='disabled')
        d_button_window = self.canvas.create_window((canvas_width-150),20,window=self.delete_button)
##Create a button to quit the program
        self.exit_program_button = tk.Button(self, text = 'Quit Program', command=self.stop_program)
        self.exit_program_button.configure(width=10)
        e_button_window = self.canvas.create_window(80,20,window=self.exit_program_button)
##Map out the proteins in a circle
        for i, protein in enumerate(collection):
            theta = 2*math.pi*(i)/(n)
            xpos = center[0] + canvas_width*math.cos(theta)/2.5
            ypos = center[1] + canvas_height*math.sin(theta)/2.5
            self.canvas.create_text([xpos,ypos],activefill="red",text=protein,tags='protein')
            

        self.line = None

        self.start_x = None
        self.start_y = None


##Create lines for the next couple functions
    def on_button_press(self, event):
        # save mouse drag start position
        self.start_list = self.canvas.find_overlapping(event.x-20,event.y+20,event.x+20,event.y-20)[0]
        self.start_x = event.x
        self.start_y = event.y

        # create lines
        self.line = self.canvas.create_line(self.x, self.y, 1, 1, fill="black", arrow=tk.LAST,state='disabled')

    def on_move_press(self, event):
        curX, curY = (event.x, event.y)

        # expand line as you drag the mouse
        self.canvas.coords(self.line, self.start_x, self.start_y, curX, curY)



    def on_button_release(self, event):
        self.end_x = event.x
        self.end_y = event.y
        self.end_list = self.canvas.find_overlapping(event.x-20,event.y+20,event.x+20,event.y-20)[0]
        connection_data.append([(self.start_list-4), (self.end_list-4)])
        self.delete_button.configure(state='normal')
        
        pass

    def remove_last(self):
        
        self.canvas.delete(self.line)
        del connection_data[-1]
        self.delete_button.configure(state='disabled')
        

    def quit(self):
        ##Save the connection data as network.npy (only able to be opened in python)
        
        ##The following line should be altered if more buttons are added to the canvas
#        connection_data = connection_data-1
        np.save("network",connection_data)
#        print(connection_data)
        root2.destroy()
        
    def stop_program(self):
    	root2.destroy()
    	sys.exit()


root2 = tk.Tk()
NetworkApp(root2).pack(side="top",fill="both",expand=True)
root2.mainloop()

############################### End draw_network #####################################

## Convert the network to a M initial matrix

######################### Begin convert_network.py ###################################

m_indeces = np.load("network.npy")
#m_indeces = connection_data
##Load the connections data as indeces to be used for the M matrix

##Assign a matrix of zeros.  For each index in the connections data, assign a 1 instead
m_initial = np.zeros((n,n))
for indeces in m_indeces:
    ##This may seem counterintuitive, but it's because the column protein influences the row protein in the M matrix
    row_index = indeces[1]
    column_index = indeces[0]
    m_initial[row_index,column_index] = 1
    

print m_initial
np.savetxt("Minitial.txt",m_initial,fmt="%i ", delimiter='     ')

########################### End convert_network #####################################

## Calculate the averages from the expressions data

########################## Begin avg_txt.py #########################################

#Initialize avg array
avg=[]
for i in range(0,tmpoint):
#   Load the data
    filename=f[i].replace(".fcs",".txt")
    data=np.loadtxt(open((filename),"rb"),delimiter=" ")
#   Any values in data that are less than zero get converted to zero
    data[data < 0] = 0
#Create a line of average values, append it to the avg array
    avg_line=np.mean(data,axis=0,dtype=np.float64)
    avg.append(avg_line)
avg = np.array(avg)
#Save it as avg.txt
avg_file = open(('avg.txt'),'w')
np.savetxt(avg_file,avg,fmt='%.8f')
avg_file.close()

############################ End avg_txt #######################################

avg_new = []
t1_new = []
t2_new = []
for i,boo in enumerate(checked):
    if boo==1:
        avg_new.append(avg[:,i])
        filename=f[0].replace(".fcs",".txt")
    	data=np.loadtxt(open((filename),"rb"),delimiter=" ")
    	t1_new.append(data[:,i])
    	filename=f[1].replace(".fcs",".txt")
    	data=np.loadtxt(open((filename),"rb"),delimiter=" ")
    	t2_new.append(data[:,i])

avg_new = np.array(avg_new)
avg_new = avg_new.transpose()
avg_new_file = open(('avg_new.txt'),'w')
np.savetxt(avg_new_file,avg_new,fmt='%.2f')
t1_new = np.array(t1_new)
t2_new = np.array(t2_new)
t1_new = t1_new.transpose()
t2_new = t2_new.transpose()
t1_new_filename = str(t[0]) + '_min.txt'
t2_new_filename = str(t[1]) + '_min.txt'
t1_new_file = open(t1_new_filename,'w')
t2_new_file = open(t2_new_filename,'w')
np.savetxt(t1_new_file,t1_new,fmt='%.2f')
np.savetxt(t2_new_file,t2_new,fmt='%.2f')
t1_new_file.flush()
t2_new_file.flush()
avg_new_file.flush()
t1_new_file.close()
t2_new_file.close()
avg_new_file.close()

os.system("./new_align.out")

os.rename('fort.61','avg.txt')
os.rename('fort.60','equaltime.txt')

time_slices = tmpoint-1
start_time = time_module.time()
for i in range(0,time_slices):
#	dir_name = 'timenew_'+str(i)
#	os.makedirs(dir_name)
	os.system("./Monte_Carlo.out")
	end_time = time_module.time()
	elapsed = end_time - start_time
	print 'elapsed time = ',elapsed
	os.rename('fort.41','Mmatrix_predict.dat')
	os.rename('fort.42','chi.txt')
	os.rename('fort.40','accept.txt')
	os.rename('fort.49','equil.txt')
	os.rename('fort.59','stdequil.txt')
	os.rename('fort.51','sigma.txt')

## Calculate the fluxes

############################ script for calc_flux.py #################################

fid = open("Mmatrix_predict.dat","rb")
m = np.loadtxt(fid)
avg = np.loadtxt('avg_new.txt')
dim = m.shape

current = []

##Calculate fluxes (formula taken from Sayak Mukherjee)
for i in range(0, dim[1]):
    row = []
    for j in range(0, dim[1]):
        current_val= m[i,j]*avg[0,j] - m[j,i]*avg[0,i]
        row.append(current_val)
    current.append(row)

#np.savetxt("flux_calculations.csv",current,delimiter=",")

########################### end calc_flux.py #######################################

## Create an output file for the fluxes

######################### script for create_flux_file.py ###############################

dim = m_initial.shape

column_1 = []
column_2 = []
column_3 = []
current = np.array(current)

column_1.append('Flux from:')
column_2.append('Flux value')
column_3.append('Flux to:')
##For each element of Minitial, if it exists, read in the flux from this connection and the two proetins involved
for i in range(0,dim[1]):
    for j in range(0,dim[1]):
        if (m_initial[i,j] == 1):
            column_2.append(current[i,j])
            column_1.append(proteins[j])
            column_3.append(proteins[i])

zipped = zip(column_1,column_2,column_3)
with open("final_fluxes.csv",'w') as f:
    writer = csv.writer(f,delimiter=',')
    writer.writerows(zipped)

######################## end create_flux_file.py ####################################

## Create an organized directory with any relevent results for the user

if not os.path.exists("Final Results"):
    os.makedirs("Final Results")

os.system('cp ./final_fluxes.csv ./Final\ Results/')
os.system('cp ./chi.txt ./Final\ Results/')
os.system('cp ./filename.txt ./Final\ Results/')
os.system('cp ./time.txt ./Final\ Results/')
os.system('cp ./avg_new.txt ./Final\ Results/')
os.system('cp ./protein_list.txt ./Final\ Results/')
os.system('cp ./*_min.txt ./Final\ Results/')



## Display the network with the fluxes back to the user

######################### script for show_flux.py ###################################

xpos=[]
ypos=[]
connection_data = np.load("network.npy")
connection_data = connection_data+1
line_id=[]

fluxes = current

class ShowNetwork(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.x = self.y = 0
        self.canvas = tk.Canvas(self, width=canvas_width, height=canvas_height)
        self.canvas.pack(side="top", fill="both", expand=True)
        self.canvas.bind("<Motion>",self.on_motion)

##Display the protein names in a circle again.        
        for i, protein in enumerate(collection):
            theta = 2*math.pi*(i)/(n)
            xpos.append(center[0] + canvas_width*math.cos(theta)/2.5)
            ypos.append(center[1] + canvas_height*math.sin(theta)/2.5)
            self.canvas.create_text([xpos[i],ypos[i]],text=protein)

##Create the lines connecting proteins, and associate a tag with each for the flux.            
        for line in connection_data:
            line_start = self.canvas.coords(line[0])
            line_start[0] = line_start[0] - (line_start[0]-center[0])/15
            line_start[1] = line_start[1] - (line_start[1]-center[1])/15
            line_end = self.canvas.coords(line[1])
            line_end[0] = line_end[0] - (line_end[0]-center[0])/20
            line_end[1] = line_end[1] - (line_end[1]-center[1])/20
            line = line-1
            flux = fluxes[line[1],line[0]]
            line_id.append(self.canvas.create_line(line_start,line_end,fill="black",arrow=tk.LAST,state=tk.DISABLED,tags=flux,width=2))
            
        self.v = tk.StringVar()
        self.v.set("Hover over a connection to see the flux")
        self.widget = tk.Label(self.canvas,textvariable=self.v)
        self.widget.pack()
        self.canvas.create_window(200,20,window=self.widget)

##Detect if current location overlaps with a line.  If so, show the flux of that line.
##If not, prompt user to hover over line.        
    def on_motion(self,event):
        if self.canvas.find_overlapping(event.x-10,event.y+10,event.x+10,event.y-10):
            var_exists = True
        else:
            var_exists = False
        if var_exists:
            self.line_found = self.canvas.find_overlapping(event.x-10,event.y+10,event.x+10,event.y-10)[0]
            tags = self.canvas.gettags(self.line_found)
            tags = str(tags)
            tags = tags[2:-3]
            self.v.set("Flux = " + tags)
            
        else:
            self.v.set("Hover over a connection to see the flux")


app2 = ShowNetwork()
app2.mainloop()

############################# end show_flux.py ###################################


## End program!
