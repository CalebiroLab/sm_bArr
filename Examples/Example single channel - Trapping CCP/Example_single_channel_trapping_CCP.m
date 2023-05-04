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
parameter_list_single_channel_trapping_CCP;

%% load the data
file_name='data_movie_example';
load([rawfolder,file_name,'.mat'],'X','Y','stack_CCP')


ch=1;% channel number
%% load data for state analysis
dataset.X=X(ch);% data must be contained in a cell
dataset.Y=Y(ch);% data must be contained in a cell
dataset.stack_binary_CCP=stack_CCP;

%% state analysis input
% list of channels
list_channel=[ch];
% list of interaction couples
list_interaction_couple=[];

%% run state analysis
tic
automated_list_of_states(dataset,file_name,to_be_computed,list_channel,list_interaction_couple,parameter,state_analysis_folder);
toc

