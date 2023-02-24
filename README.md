# Analysis of b2AR and bArrestin2 interactions from single-particle tracking data in CHO cells
**Plasma membrane preassociation drives b-arrestin coupling to receptors and activation**  
Jak Grimes, Zsombor Koszegi, Yann Lanoiselée, Tamara Miljus, Shannon L. O’Brien, Tomasz M. Stepniewski, Brian Medel-Lacruz, Mithu Baidya, Maria Makarova, Ravi Mistry, Joëlle Goulding, Julia Drube, Carsten Hoffmann, Dylan M. Owen, Arun K. Shukla, Jana Selent, Stephen J. Hill, Davide Calebiro  
bioRxiv, 2022.11. 15.516577

## System requirements
MATLAB R2018b
## Analysis requirements
#### Movie Names
Movies of Receptor and Arrestin are recorded simultaneously as .tif files. Receptor and arrestin movies should have the same base name (example ‘{movies_name_basis}1’) but with a suffix ‘-C1’ for receptor and a suffix ‘-C2’ for arrestin.
#### Folders
The file ‘global_folders.m’ contains all the paths to the folders used in the analysis. The user should create all the folders specified in this file and update the paths in the file accordingly. All raw data should be stored in the folder specified as 'global_folders.rawfolder'.
#### Movie list
A file should be created that contains all movies organized by groups inside a MATLAB structure as follows:
movie_list.GROUP_NAME_1=['movie1 movie2 movie3'];
movie_list.GROUP_NAME_2=['movie4 movie5 movie6'];

To these movie names, suffixes ‘-C1’ and ‘-C2’ will be attached automatically for receptor and arrestin respectively during analysis.

####	Analysis parameters
The structure parameter contains all the parameters used for analysis. It can be loaded by calling the script ‘parameter_list.m’

#### Coordinate alignment
A .mat file should be created containing 2 variables of 1x3 cells containing ‘tform’ transformations for channel alignement. The first variable named ‘t_piecewise_linear’ should contain nothing in first cell, the tform for transformation from C2 to C1 in second cell and the tform for transformation from C3 to C1 in the third cell. The second variable ‘t_piecewise_linear_rev’ has the similar structure but with reverse transformations (empty, C1 to C2, C1 to C3).  
We recommend using a 'piecewise linear' tform although an 'affine' tform transformation could be used as well. 
The alignment matrix should be stored in the folder specified as ‘global_folders.rawfolder’.
## Trajectory analysis
#### Detection and tracking
Detection and tracking are done using the u-track software, which can be downloaded here:
https://github.com/DanuserLab/u-track
>Jaqaman, K., Loerke, D., Mettlen, M. et al. Robust single-particle tracking in live-cell time-lapse sequences. Nat Methods 5, 695–702 (2008). https://doi.org/10.1038/nmeth.1237

The u-track parameters used in our study for detection and tracking can be found in the script called ‘detection_tracking_parameters.m’. 
After excluding trajectory coordinates that lie outside the region of interest, the output of u-track is converted to a custom format using the script ‘utrack_to_OBJ_v1.m’, which reorganises compound trajectories from u-track into separate objects.  
The script is adapted from 
>Cocucci E. et al. The First Five Seconds in the Life of a Clathrin-Coated Pit. Cell 150, 3 495-507 (2012).
This will generate an intermediate file ‘{movie_name_basis}-C{n}_gui2.mat’ a the file containing trajectories in our format named ‘{movie_name_basis}-C{n}_gui2_steps.mat’ in the folder ‘global_folders.rawfolder’.  
‘{movie_name_basis}-C{n}_gui2.mat’ file contains three variables:
1. OBJ: a structure where data coordinates are contained in OBJ.xR and OBJ.yR and spot intensities in OBJ.trjRB, when a molecules disapears the last coordinate is repeated for all subsequent frames.
2. OBJ2: a structure where data coordinates are contained in OBJ2.xR and OBJ2.yR and spot intensities in OBJ2.trjRB, when a molecules disapears all subsequent frames are NaN.
3. IFO: structure containing informations on the data.
4. TME: variable containing the time point corresponding to each frame

This file contains several matrices. For each movie, a new field ‘IFO.calmatrix’ should be  appended to the IFO structure, which contains the file name of the alignment matrix (Example: IFO.calmatrix=’alignement_matrix’).

#### Interaction analysis
The interaction analysis is performed using the script ‘cycles_interaction.m’. The script computes the colocalization events between  C1 and C2 as described in Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

This script generates the following files in the folder ‘global_folders.rawfolder’:
1. ‘{movie_name_basis}-C1{movie_name_basis}-C2_intmatrix2.mat’: contains the 'interaction matrix', a n_i x 10 matrix classifying all interations found, where First column and Second column are the first and last frame of interaction, column 3 and 7 are the channel numbers, column 4 and 8 correspond to the trajectory number in OBJ2 of the GUI2_steps file.   

2. ‘{movie_name_basis}-C1{movie_name_basis}-C2_intmatrix_0pShiftX-0pShiftY-0fShiftT.mat’: contains a structure 'Ch' that contains all the trajectories post alignement combined with the interaction matrix and a list of all interaction start, and a list of interaction 'families', portions of interactions that together form a longer interaction.

The algorithm first does gap closing of the X and Y coordinates, then it computes the interaction matrix. Then an optimisation routine is used to link interactions portions together into longer 'interaction families' by trying to minimize the number of interaction families and starting from the case where each interaction is in a different family. This process is necessary to deal with situations with merges and splits. For example, if two receptor molecules are at the same place, their trajectory will exactly the same until they split. If while merged the receptors colocalise with an arrestin molecule then both receptors with be recorded to colocalise with the same arrestin in the 'interaction matrix'. Because an arrestin cannot interact with two receptors at the same time, the optimisation routine finds to which receptor the interaction should be attributed depending on which continues to colocalise after the split.
Two interactions can be connected only if both interactions contain the same molecules 



Please cite 
>Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).




#### Analysis of diffusive states
The trajectories are first analysed with an algorithm that detects transient trapping events for each trajectory and each channel, described in 
>Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021). 

The code for transient trapping detection can be downloaded here: 
https://github.com/YannLanoiselee/Transient_trapping_analysis

Please cite 
>Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021).

Then, the information about transient trapping is combined with that of colocalization between C1 and C2 over time using the script “cycle_states_forced_or_not.m”. This generates a file ‘{movie_name_basis}-C{n}_list_states.mat’ for each movie and channel n in the folder ‘global_folders.state_analysis_folder’.

#### Time-averaged MSD
The TAMSD is computed for each trajectory in C1 and C2 using the function ‘cycle_TAMSD_GABAB_Filamin.m’. For each trajectory, the analysis estimates the anomalous exponent α and the generalized diffusion coefficient D_α by fitting the TAMSD curve as a function of lag-time with the formula for the average TAMSD in 2 dimensions:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\langle%20\delta^2(\Delta,t)\rangle=4D_\alpha%20\Delta^\alpha+4\sigma^2"/>

The function generates a result matrix saved in ‘{movie_name_basis}-C{n}_MSD.mat’ containing the computed generalized diffusion coefficient and the anomalous exponent values in the first and second column, respectively.


