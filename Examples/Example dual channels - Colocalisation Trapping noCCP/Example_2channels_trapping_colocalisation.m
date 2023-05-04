% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
clear;
close all;

%% define load and save paths
% define load path as the folder where the script is located 
script_path = mfilename('fullpath');
idx=find(ismember(script_path,filesep),1,'last');
rawfolder=script_path(1:idx);
% folder in which interaction analysis data is saved
interaction_data_folder=rawfolder;
% folder in which list of states is saved
state_analysis_folder=rawfolder;

%% load the script that contains all analysis parameters and states definition
parameter_list_2channels_trapping_colocalisation;

%% load the data
file_name='data_movie_example';
load([rawfolder,file_name,'.mat'],'X','Y') % data containing trajectory cooridnates for the two channels

%% Gap closing
gap_method={'linear_interpolation','brownian_bridge'};
n_method=1;
disp(['Gap closing'])
tic
[X,Y]=gap_close(X,Y,gap_method{n_method});
toc

%% Interaction analysis
disp(['Computing interaction matrix'])
tic
[int_matrix]=interaction_matrix(X,Y, 1,parameter.FrameNumber,parameter.searchRadius_pixel);
toc
disp(['Maximizing interactions'])
tic
gap_wind=3;
[interaction_data] = compound_interactions_maximize(int_matrix,gap_wind);
toc
save([interaction_data_folder,filesep,file_name,'_interaction_data.mat'],'interaction_data','int_matrix');
disp(['Compute swapping'])
tic
[X,Y,int_matrix,interaction_data] = XY_swap_after_interaction(X,Y,int_matrix,interaction_data);
toc
interaction_data.int_matrix=int_matrix;
save([interaction_data_folder,filesep,file_name,'_swapped.mat'],'X','Y','interaction_data','int_matrix');

%% load data for state analysis
dataset.X=X;
dataset.Y=Y;

dataset.interaction_data=interaction_data;
dataset.interaction_data.int_matrix=int_matrix;

%% state analysis input
% list of channels
list_channel=[1,2];
% list of interaction couples
list_interaction_couple=[1,2];

%% run state analysis
tic
automated_list_of_states(dataset,file_name,to_be_computed,list_channel,list_interaction_couple,parameter,state_analysis_folder);
toc

