# Xsorb

xSorb is a lightweight, flexible and user-friendly tool to generate and run ab initio calculations of molecules adsorbed over surfaces. 
Relying on either Quantum Espresso or VASP for ab initio calculations and on two well-established Python libraries (ASE and Pymatgen) to generate the structures, it allows to produce and explore a large number of adsorption configurations of a molecule on a surface, finally identifying the most stable configuration  and computing the associated adsorption energy.

# Installation

## Download and setup

The installation procedure for Xsorb is easy and straightforward since it is Python-based. It is first necessary to download the program by cloning this repository into the desired local machine:  

`git clone https://gitlab.com/triboteam/xsorbed.git`

Once the download is completed, it is necessary to go to iside the downloaded folder and run: 

`bash install.sh`

to add the executable in the `PATH` variable. With this operation, the user can launch the program from any folder. The next step is to install the required dependencies stored in the 'requirements.txt' file. This step can be performed with the command:

`python3 -m pip install -r requirements.txt`

Xsorb is interfaced with Quantum Espresso and VASP to perform the DFT calculations. It is, therefore, necessary to compile and install the corresponding executable before launching the calculations.

After the installation, the program works through a command-line interface (CLI) by running:  


`xsorb [command] [parameters]`


Where `command` specifies which operation to perform and `parameters` are additional parameters available for specific commands.


## Dependencies

The program requires the following Python libraries:
- [`ase`](https://wiki.fysik.dtu.dk/ase/index.html)
- [`pymatgen`](https://pymatgen.org/)
- [`numpy`](https://numpy.org/)
- [`matplotlib`](https://matplotlib.org/)
- [`pandas`](https://pandas.pydata.org/)


### Optional requirements
In order to generate very high quality images and animations of the atoms inside the simulation cell, the program relies on the external program [`povray`](http://www.povray.org/), a ray-tracing based renderer, which on Linux systems can be installed directly from command line (e.g. using apt: `sudo apt install povray`). 
 
For the animations with povray, the [`ffmpeg`](https://ffmpeg.org/) package is also needed, which again can be simply installed on Linux system with `sudo apt install ffmpeg`

Note however that this is not mandatory, and you can still generate static images and animations by using ase's own renderer, although with a lower image quality.

# Program use

The program works through a command-line interface (CLI) by running:  

`xsorb [command] [parameters]`

where `command` specifies which operation to perform and `parameters` are additional parameters available for specific commands.

A comprehensive list of all the commands with a brief description is present in the program's help, accessible by typing:

`xsorb -h` or `xsorb --help`

The program requires two input files containing the (optimized) structures of the molecule and the slab that the user wants to study. All the options needed to generate the adsorption configurations and the computational parameters must be included in a 'settings.in' file, a template of which is provided in the main folder of this repository. Two versions are present: one with all the parameters, each with a brief description, and another, called 'settings_minimal.in', with only the mandatory options.

## Main workflow

The usage procedure of the program is constituted by three main steps:

i) automatic generation of many adsorption configurations through the molecule rotation and translation with respect to the surface; 

ii) identification of the most relevant ones by screening their energies through a partial optimization; 

iii) full structural optimisation of the configurations identified as most relevant. 

A detailed description of the worfklow and two examples of usage can be found in the repository's WIKI (https://gitlab.com/triboteam/xsorbed/-/wikis/home) and in the article related to this program (https://doi.org/10.1016/j.cpc.2023.108827).
