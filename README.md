# MESA_from_TAPIR
Fortran Subroutine to use a standalone installation of MESA inside the TAPIR code.

Before this subroutine is able to work, you have to install [MESA](http://mesa.sourceforge.net/prereqs.html) and should test it.

# How this works
The MESA working directory is set up to run a piecewise MESA calculation. This is steered by the programm `<hello>` in *MESA_from_TAPIR.f90* by the use of the `<start_MESA>` subroutine.
The subroutine is fed the starting and end point of the calculation as well as the start and end accretion rate. It will then call *rn_tapir*, which starts the MESA run in one of two ways:
* If the start time is 0, it will call `<rn>` which reads in the initial model and evolves it until `<time_end>` with an accretion rate that is linearly interpolated between `<rate_start>` and `<rate_end>`.
Once it has reached `<time_end>`, MESA will save a binary file (also called photo in MESA) from which we can later resume the calculation.
* If the start time is > 0 it will call `<re>` which starts from the last photo and evolves it until `<time_end>` with an accretion rate that is linearly interpolated between `<rate_start>` and `<rate_end>`.
Once it has reached `<time_end>. MESA will save a binary file (also called photo in MESA) from which we can later resume the calculation.

At every timestep MESA writes out information to *LOGS/history.data* that include e.g. effective temperature, luminosity, star_mass, etc. To talk with the *MESA_from_TAPIR.f90* subroutine, MESA will write the 
wanted information to *MESA_for_TAPIR.txt* in each timestep. In our case this includes the effective temperature and the time. After the MESA run has ended, `<start_MESA>` will read in the the informatin from this file.
If the time does not match `<time_end>` and error is raised.

For the `<rn_tapir>` file to work on your machine you will need to change the environmental variables set in the first few lines accordingly. This is how it is setup for my machine.:

```
export MESA_DIR=~/MESA/mesa-r12778 (Path to MESA installation)
export OMP_NUM_THREADS=16 (Amount of course to run MESA on)
export MESASDK_ROOT=~/mesasdk (Path to mesasdk installation)
```


This repository includes a MESA working directory and the subroutine to include MESA in TAPIR. Here is a breakdown of what files are doing what:


* **The MESA work directory**:
  * *inlist*: is the basic inlist that MESA star will read at the initialisation. It points towards inlist_project, the file with all controls
  * *inlist_project*: includes the controls for the MESA evolution
  * *inlist_pgstar*: includes the controls for th PGPLOT window if run. Including this is set to off.
  * *history_columns.list*: defines the output of the MESA run in LOGS/history.data, a file that stores overall properties of every model timestep such as log_Teff, etc
  * *initial_model_h2_birthline*: The initial model for this trial run. A 1 Msun model that has been created using create_pre_main_sequence_model and evolved until the central H2 is below 1d-10
  * *mk*: bash file to compile this MESA work directory to create the Star executable
  * *rn*: bash file to run MESA star in the current directory
  * *re*: bash file to resume MESA run in the current directory
  * *clean*: bash file to clean up things before compiling
  * **src**:
    * *run_star_extras.f90*: This file allows for the modificatio of MESA. It included for example the treatmen of accretion.
  * **make**: folder that includes stuff for compiling
 
* *accr_rate_burst.dat*: File that contains accretion rates
* *accr_rate_noburst.dat*: File that contains accretion rates
* *MESA_for_TAPIR.txt*: text file over which MESA talks to the subroutine. Includes the data that is fed back to TAPIR. In this stage the effective temperature and final time (for error control).
* *MESA_from_TAPIR.f90*: Subroutine that runs MESA in TAPIR
* *output.txt*: text file that hold the command line output of the last MESA run. For debugging.
* *restart_photo*: some file that MESA uses for running restarts. not sure if we need this.
* *rn_tapir*: bash file that runs MESA evolution from time_start to time_end with accretion rates linearly interpolated from rate_start to rate_end
* *inlist_project_file*: blueprint for inlist_project. Is copied to inlist_project from rn_tapir before some values are replaced.
