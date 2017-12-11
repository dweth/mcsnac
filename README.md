# single-cell-snac

The Single Cell Signaling Network Analysis Code, or Single Cell SNAC, is a set of files designed to allow anyone to analyze signaling network kinetics of interest from large mass cytometry datasets.  It utilizes Simulated Annealing to propose a first-order reaction network between molecules of interest.  NOTE: THIS CODE IS ONLY CURRENTLY FUNCTIONAL ON MAC OS X.

# Table of Contents
*[Quick Start](#quick-start)
*[Repository Contents](#repository-contents)
*[Installing Gfortran](#installing-gfortran)
*[Detailed Tutorial](#detailed-tutorial)

# Quick Start
1. Make sure you have gfortran downloaded (see "Installing Gfortran").
2. Open a Terminal (/Applications/Utilities/Terminal).
3. Drag and drop the downloaded file "drag_me_to_terminal.py" to the Terminal and press Enter.
4. Follow the program's instructions.

# Repository Contents (Eventually Needs Edit)
The documents contained in this repository include:
1. Before Running
   -These are instructions for installing and downloading all the necessary prerequisites to run the software.
2. Instructions
   -These are instructions as to how to use the code once you have the necessary prerequisites.
3. drag_me_to_terminal.py
   -This file is the main file of the program.  It executes much of the code, and calls the other scripts in the repository.
4. new_align.f
   -This file organizes the necessary files in a new directory.  It "preps" for the Simulated Annealing scheme.
5. Monte_Carlo.f
   -This file executes the Simulated Annealing scheme.
   
# Installing Gfortran
Install gfortran using the following link: https://gcc.gnu.org/wiki/GFortranBinaries

You’ll see several options for installers.  Select the one that corresponds to your operating system and follow all instructions for installation.

INSERT LINKS TO PHOTOS HERE.

# Detailed Tutorial
*[Before Starting](#before-starting)
*[Running the Program](#running-the-program)

# Before Starting
Make sure all of the following are true before you utilize the software:

1. You have downloaded gfortran.  See "Installing Gfortran".
2. You have 2 .fcs files that you would like to analyze.
   2b. You know the time after t=0 the data were collected (e.g. 8 minutes and 16 minutes).
3. You have a signaling network to analyze in mind.
4. You have time.  The analysis can take anywhere from a couple hours to several days.  You can still use the computer during this time for other tasks, but keep it plugged in to a power source.

Now you are ready to start the program!  The instructions below detail how to do this.

# Running the Program

1.	Open a Terminal.  This is done by opening the Finder, going to the Applications folder, opening the Utilities folder, and double clicking on “Terminal”. 

A screen similar to that shown below should appear.  This is the Terminal.

2.	Locate the file titled “drag_me_to_terminal.py” that was downloaded from GitHub.  Drag and drop this file into the Terminal.

3.	Go to the Terminal and press the “Enter” key.  This will begin running the software. (From here on, the software will guide you through using it)

4.	Provide the prompt with the two times that data was collected at.

5.	Select each of the data files to analyze.

6.	Select the proteins to construct the network from.

7.	Draw the network.

8.	The simulation will begin.  To see the status of it, the bottom line displayed in the Terminal should show “Run [X] out of 2000”.  The simulation will complete when it reaches 2000.

9.	Once complete, you may view the values of the flux between network connections by hovering over the arrows you drew previously.

10.	The final results are created in a folder named after the time you began the simulation.  The flux values can be found in the file titled “final_fluxes.csv”.  



If desired, these can be visualized fairly easily using the free, open-source software CytoScape. 
