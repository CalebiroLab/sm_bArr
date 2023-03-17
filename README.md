# Single molecule analysis of protein-protein spatiotemproal co-dynamics
**Plasma membrane preassociation drives b-arrestin coupling to receptors and activation**  
Jak Grimes, Zsombor Koszegi, Yann Lanoiselée, Tamara Miljus, Shannon L. O’Brien, Tomasz M. Stepniewski, Brian Medel-Lacruz, Mithu Baidya, Maria Makarova, Ravi Mistry, Joëlle Goulding, Julia Drube, Carsten Hoffmann, Dylan M. Owen, Arun K. Shukla, Jana Selent, Stephen J. Hill, Davide Calebiro  
bioRxiv, 2022.11. 15.516577

[![DOI](https://zenodo.org/badge/570954506.svg)](https://zenodo.org/badge/latestdoi/570954506)

## System requirements
MATLAB R2018b  
Image Processing Toolbox  
Statistics and Machine Learning Toolbox

## Analysis requirements


## Steps to follow for analysis


#### Clathrin Coated-Pits binary mask
CCP movies are not tracked because CCP are larger than diffraction limit. Instead a binary mask is made using the function 'binary_msk_CCP'.

#### Interaction analysis
The interaction analysis is performed using the script ‘cycles_interaction.m’. The script computes the colocalization events between C1 and C2 as described in Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

A file is created  where the trajectory are reorganised according to the interaction linking performed by the optimisation algorithm that also contains the interaction information.

Please cite 
>Sungkaworn, T. et al. Single-molecule imaging reveals receptor-G protein interactions at cell surface hot spots. Nature 550, 543–547 (2017).

#### Analysis of diffusive states
The trajectories are first analysed with an algorithm that detects transient trapping events for each trajectory and each channel, described in Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021). 

The code for transient trapping detection can be downloaded here: 
https://github.com/YannLanoiselee/Transient_trapping_analysis

Please cite 
>Lanoiselée, Y., Grimes, J., Koszegi, Z. & Calebiro, D. Detecting transient trapping from a single trajectory: A structural approach. Entropy 23, 1–16 (2021).

Then, the information about transient trapping is combined with that of colocalization between recpetor and arrestin as well as with the information about colocalisation of receptor or arrestin with Clathrin Coated Pits over time using the script “cycle_states_forced_or_not.m”. This generates a file ‘{movie_name_basis}-C{n}_list_states.mat’ for each movie and channel n in the folder ‘global_folders.state_analysis_folder’.
<!---
#### Time-averaged MSD
The TAMSD is computed for each trajectory in C1 and C2 using the function ‘cycle_TAMSD.m’. For each trajectory, the analysis estimates the anomalous exponent α and the generalized diffusion coefficient D_α by fitting the TAMSD curve as a function of lag-time with the formula for the average TAMSD in 2 dimensions:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\langle%20\delta^2(\Delta,t)\rangle=4D_\alpha%20\Delta^\alpha+4\sigma^2"/>

The function generates a result matrix saved in ‘{movie_name_basis}-C{n}_MSD.mat’ containing the computed generalized diffusion coefficient and the anomalous exponent values in the first and second column, respectively.
--->

