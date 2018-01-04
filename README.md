# single-cell-snac

The Single Cell Signaling Network Analysis Code, or Single Cell SNAC, is a set of files designed to allow anyone to analyze signaling network kinetics of interest from large mass cytometry datasets.  It utilizes Simulated Annealing to propose a first-order reaction network between molecules of interest.  NOTE: THIS CODE IS ONLY CURRENTLY FUNCTIONAL ON MAC OS X.

## Table of Contents
* [Quick Start](#quick-start) 
* [Repository Contents](#repository-contents) 
* [Installing Gfortran](#installing-gfortran) 
* [Detailed Tutorial](#detailed-tutorial) 
* [Viewing Results After Completion](#viewing-results-after-completion) 

## Quick Start
1. Make sure you have gfortran downloaded (see "Installing Gfortran").
2. Open a Terminal (/Applications/Utilities/Terminal).
3. Drag and drop the downloaded file "drag_me_to_terminal.py" to the Terminal and press Enter.
4. Follow the program's instructions.

## Repository Contents
The documents contained in this repository include:
1. README
   -Contains instructions for downloading necessary prerequisites, running the software, and other notes.
2. drag_me_to_terminal.py
   -This file is the main file of the program.  It executes much of the code, and calls the other scripts in the repository.
3. new_align.f
   -This file organizes the necessary files in a new directory.  It "preps" for the Simulated Annealing scheme.
4. Monte_Carlo.f
   -This file executes the Simulated Annealing scheme.
   
## Installing Gfortran
Install gfortran using the following link: https://gcc.gnu.org/wiki/GFortranBinaries

You’ll see several options for installers.  Select the one that corresponds to your operating system.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/gfortran_binaries.png)

Once gfortran is downloaded, open the .dmg file from your downloads.  This will contain a file 'gfortran.pkg'.  Open this as well.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Downloads_screen.png)
![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Gfortran_pkg.png)

This will take you to instructions for installation.  Accept any terms and install in the suggested location.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Gfortran_Instruct.png)

Once you see this screen, gfortran has been installed and you are ready to move on!

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Gfortran_Installed.png)


## Detailed Tutorial
* [Before Starting](#before-starting) 
* [Running the Program](#running-the-program)

### Before Starting
Make sure all of the following are true before you utilize the software:

1. You have downloaded gfortran.  See [Installing Gfortran](#installing-gfortran).
2. You have 2 .fcs files that you would like to analyze.
3. You know the time after t=0 the data were collected (e.g. 8 minutes and 16 minutes).
4. You have a signaling network to analyze in mind.
5. You have time.  The analysis can take anywhere from a couple hours to several days.  You can still use the computer during this time for other tasks, but keep it plugged in to a power source.

Now you are ready to start the program!  The instructions below detail how to do this.

### Running the Program

1.	Open a Terminal.  This is done by opening the Finder, going to the Applications folder, opening the Utilities folder, and double clicking on “Terminal”.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Open_Terminal_1.png)
![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Open_Terminal_2.png)
![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Open_Terminal_3.png)

A screen similar to that shown below should appear.  This is the Terminal.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Empty_terminal.png)

2.	Locate the file titled “drag_me_to_terminal.py” that was downloaded from GitHub.  Drag and drop this file into the Terminal.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/drag_me_to_terminal.png)
![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Terminal_with_path.png)

3.	Go to the Terminal and press the “Enter” key.  This will begin running the software. (From here on, the software will guide you through using it)

If you see the following screen, click the "Install" button to install Command Line Developer Tools.  This may happen the first time you run the program; if so, simply re-run it by pressing the [up] key and then [enter] in the Terminal.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Developer_Tools_Install.png)

4.	Provide the prompt with the two times that data was collected at (the example collected data at 8 and 16 minutes).

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_1b.png)

5.	Select each of the data files to analyze.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_1c.png)
![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_2.png)

6.	Select the molecules to construct the network from by clicking the checkboxes next to the molecules of interest.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_6.png)

7.	Draw the network by clicking and holding on one molecule, then dragging an arrow to the molecule to connect it to.  Please note the direction of the arrow.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_9.png)

When you're done, click the "Finish" button.

8.	The simulation will begin.  To see the status of it, the bottom line displayed in the Terminal should show “Run [X] out of 2000”.  The simulation will complete when it reaches 2000.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_10.png)

9.	Once complete, you may view the values of the flux between network connections by hovering over the arrows you drew previously.

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_12.png)

10. The final results are created in a folder named after the time you began the simulation.  The flux values can be found in the file titled “final_fluxes.csv”.


## Viewing Results After Completion

When you run the software, it automatically creates a folder of the exact datetime the program began.  For example, if the program was started at 1:05:30 p.m. on October 20th, 2016, you would see a folder titled something like "10_20_2016_13_05_30" in the folder containing the file drag_me_to_terminal.py.  If you enter this folder, you will find another folder titled "Final Results".  Enter this to find the following documents.
* [Fluxes](#fluxes)
* [Chi Squared](#chi-squared)
* [Filenames Used](#filenames)
* [Times Used](#times-used)
* [Data](#data)

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_13.png)

### Fluxes
If you enter this folder, you will find a file "final_fluxes.csv" for any completed computation.  This contains the results of the analysis in a format readable by Microsoft Excel.  If desired, the results in "final_fluxes.csv" can be visualized fairly easily using the free, open-source software [CytoScape].

![alt tag](http://planetx.nationwidechildrens.org/~jayajit/software_pics/Tutorial_15.png)

### Chi Squared
To see the progression of the chi-squared value as the simulation executed, open the file titled "chi.txt".  The final value in this file is the chi squared associated with the final results.

### Filenames
If you forget which data files you analyzed, the filenames and their paths can be found in the file "filename.txt".  
### Times Used
Additionally, the times entered by the user for the data can be found in "time.txt".  

### Data
The data itself for the proteins of interest can be viewed by opening the files named after the times (e.g. 8_min.txt and 16_min.txt).  The averages of these data can be seen in the file "avg_new.txt".  The order of these data columns is in the order they appear in the file "protein_list.txt".

[CytoScape]: http://www.cytoscape.org/
