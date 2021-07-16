*****************************************************************************************************************
***************************************************** ABOUT *****************************************************
*****************************************************************************************************************
The sibulator.py script takes in several command-line arguments for different details of its sibling simulation.
It can be run in several different ways and outputs a file with the specified name and type (type is given by
the extension type i.e. .csv, .txt, etc.). For best readability, write to a .csv file. The first line of the
output will be the known contributor profile and all subsequent ones are randomly generated.

The sibling profiles are generated as follows:
Input: known sample, allele frequency data, subpopulation, N (number of random profiles), file to write to
For each n = 1,...,N:
	- Create an empty profile
	For each location listed in the known profile:
		- Randomly choose to replace 0, 1, 2 of the alleles at the location with probabilies 0.25, 0.5, 0.25,
		  respectively
			- If replacing only one, randomly choose which one of the alleles to replace
		- Query the allele frequency database specified (nist or strider) for the allele frequencies at given
		  location
		- Sample the requisite number of alleles given the allele frequencies for the current location
		- Add the location with the chosen alleles to the current profile
	End for
End for
Return N generated samples and write them to the specified file

*****************************************************************************************************************
*************************************************** HOW TO RUN **************************************************
*****************************************************************************************************************

There are a few entry points into the sibulator program

-------------------------------------------- Direct via Command Line --------------------------------------------
Script sibulator.py takes 5 following command-line arguments for parameters on how to generate sibling samples:
	1) known_sample_path - local path to file of known sample
	2) allele_freq_data - which allele frequency to use (nist or strider)
	3) subpop - which subpopulation's allele frequencies to use (see NIST 
				and STRidER documentation for available subpopulations)
	4) num_siblings - number of sibling samples to generate
	5) dest_path - name of file to create and store resulting samples

Below is example line to run in command line/terminal (replace non-"--" values for own use)

python sibulator.py --known_sample_path AfAm\ Noncontributor\ Database.csv --allele_freq_data nist --subpop AfAm --num_siblings 1000 --dest_path test.csv

*NOTE* This requires your working directory of the terminal/command line prompt window to be the location of the
code i.e. must be pointing to the src folder. Take the following steps to accomplish this:
	1) Open a terminal/command-line prompt window 
		- Mac: command+space, type Terminal into search bar and press enter
		- Windows: search for command prompt in start menu
	2) Open a finder/files window and go to where you've downloaded the "src" folder
	3) Type "cd " and click and drag the "src" folder into the terminal/command line prompt window
		- This should autofill the path to the downloaded folder
	4) Hit enter
	5) Now you can run the command line argument given above (also use this process for the GUI via command line)

---------------------------------------------- GUI via Command Line ----------------------------------------------

Can also be used via the GUI script sibulator_gui.py. Note this is primitive and all values must be exactly correct (path names, subpopulation name are case sensitive) on input will have to restart window

GUI started using the following command-line entry (must be in the local path as described in NOTE above)

python sibulator_gui.py

---------------------------------------------- GUI via Double Click ----------------------------------------------

On Windows (untested): double-click Sibulator.bat

On MacOS (confirmed works): double-click Sibulator