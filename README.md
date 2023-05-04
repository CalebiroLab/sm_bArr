# Single molecule analysis of protein-protein interactions and spatiotemporal dynamics
This repository contains the MATLAB scripts used for the single-moelcule analyses on receptors and beta-arrestins in:
**Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation**  
Jak Grimes, Zsombor Koszegi, Yann Lanoiselée, Tamara Miljus, Shannon L. O’Brien, Tomasz M. Stepniewski, Brian Medel-Lacruz, Mithu Baidya, Maria Makarova, Ravi Mistry, Joëlle Goulding, Julia Drube, Carsten Hoffmann, Dylan M. Owen, Arun K. Shukla, Jana Selent, Stephen J. Hill, Davide Calebiro  
Grimes, Koszegi, Lanoiselée et al., 2023, Cell 186, 1–18
https://doi.org/10.1016/j.cell.2023.04.018

[![DOI](https://zenodo.org/badge/570954506.svg)](https://zenodo.org/badge/latestdoi/570954506)

## System requirements
MATLAB R2018b  
Image Processing Toolbox  
Statistics and Machine Learning Toolbox

## Analysis Inputs
Trajectories 'X' and 'Y' coordinates must be stored in two MxN matrices where 'M' is the number of trajectories and 'N' is the number of frames. Trajectory coordinates must be in pixel units. Missing data points must be 'NaN' values.
Clathrin-Coated Pit (CCP) localizations must be stored either as a binary stack of size LxLxN where 'L' is the number of pixels (assuming square images) or as a LxL single image, where 1 denotes the presence of CCP and 0 denotes its absence.

The CCP pixel size must match that of the trajectory coordinates.

## Steps to follow for the analysis

#### Interaction analysis
The interaction analysis computes the colocalization events between channel 1 (C1)  and channel 2 (C2) as described in Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

An output file is created that contains the results of the interaction analysis and where the trajectories are reorganized to optimize the results by linking interaction fragments that are contiguous in space and time. 

Please cite 
>Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

#### Analysis of diffusive states
The trajectories are first analysed with an algorithm that detects transient trapping events for each trajectory and each channel, described in Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021). 

The code for transient trapping detection can be downloaded here: 
https://github.com/YannLanoiselee/Transient_trapping_analysis

Please cite 
>Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021).

Then, the information about transient trapping is combined with that about colocalization between each molecule in C1 and C2 and their localization in CCPs.

The analysis generates an output file ‘{movie_name_basis}-C{n}_list_states.mat’ for each channel {n}.

## State assignment 
The analysis generates separate binary matrices that contain information about transient trapping, interaction with molecules in the other channel and localization in CCPs for each molecule at each frame. This information is used to build a matrix that contains the state assign to each moleucle at each frame. A dull state is added to represent molecules before/after their appearance/disappearance, for instance due to movement between the cytoplasm and the plasma membrane. 

Depending on the settings, the scripts can be used to obtain different type of information.

## Examples

#### One channel only with CCPs
In this example, the script computes the transient trapping of trajectories, the presence/absence of CCP localization and combines this information in 3+1 states stored in a MxN matrix. 

#### Two channels

In this example, the script starts with the colocalization between the molecules in the two channels, then computes their transient trapping and combines the information about the resulting 4+1 states in a MxN matrix. 

#### Two channels with CCPs

In this example, the script starts with the colocalization between the molecules in the two channels, then computes their transient trapping and presence/absence of CCP localization and combines the information about the resulting 6+1 states in a MxN matrix. 
