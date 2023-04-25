# Single molecule analysis of protein-protein spatiotemproal co-dynamics
**Plasma membrane preassociation drives b-arrestin coupling to receptors and activation**  
Jak Grimes, Zsombor Koszegi, Yann Lanoiselée, Tamara Miljus, Shannon L. O’Brien, Tomasz M. Stepniewski, Brian Medel-Lacruz, Mithu Baidya, Maria Makarova, Ravi Mistry, Joëlle Goulding, Julia Drube, Carsten Hoffmann, Dylan M. Owen, Arun K. Shukla, Jana Selent, Stephen J. Hill, Davide Calebiro  
bioRxiv, 2022.11. 15.516577

[![DOI](https://zenodo.org/badge/570954506.svg)](https://zenodo.org/badge/latestdoi/570954506)

## System requirements
MATLAB R2018b  
Image Processing Toolbox  
Statistics and Machine Learning Toolbox

## Analysis Inputs
Trajectories 'X' and 'Y' coordinates must be stored in two MxN matrices where 'M' is the number of trajectories and 'N' is the number of frames. Trajectory coordinates must be in pixel units. Missing data points must be 'nan' values.
Clathrin-Coated Pit (CCP) localisations must be stored either as a binary stack of size LxLxN where 'L' is the number of pixels (assuming square images) or as a LxL single image, where 1 denotes the presence of CCP and 0 denotes its absence.

CCP pixel size must match trajectories coordinates pixel size.

## Steps to follow for analysis

#### Interaction analysis
The interaction analysis computes the colocalization events between C1 and C2 as described in Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

A file is created  where the trajectory are reorganised according to the interaction linking performed by the optimisation algorithm that also contains the interaction information.

Please cite 
>Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

#### Analysis of diffusive states
The trajectories are first analysed with an algorithm that detects transient trapping events for each trajectory and each channel, described in Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021). 

The code for transient trapping detection can be downloaded here: 
https://github.com/YannLanoiselee/Transient_trapping_analysis

Please cite 
>Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021).

Then, the information about transient trapping is combined with that of colocalization between molecule 1 and and arrestin as well as with the information about colocalisation of receptor or arrestin with Clathrin Coated Pits over time. This generates a file ‘{movie_name_basis}-C{n}_list_states.mat’ for each channel {n}.

## States definition 
The analysis generates independent binary vectors for transient trapping, interaction with another molecule type,   
## Examples

#### 1 molecule type with CCP
In this case the analysis computes the transient trapping of trajectories, the presence or not of molecules at CCP and combines this information in 3+1 states stored in a MxN matrix. 

#### 2 molecule types

In this case the analysis start with the colocalisation between the two molecules, then computes the transient trapping for trajectories of each molecule type and combines this information in stored in a MxN matrix. 

#### 2 molecule types with CCP

In this case the analysis start with the colocalisation between the two molecules, then computes the transient trapping for trajectories of each molecule type, computes the presence or not of molecules at CCP and combines this information into 6+1 states. 

