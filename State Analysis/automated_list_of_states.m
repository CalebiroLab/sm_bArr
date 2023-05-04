% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_not_working]=automated_list_of_states(dataset,moviename,to_be_computed,list_channel,list_interaction_couple,parameter,save_path)

% dataset : structure that contains informations for interaction analysis 
% dataset.X: trajectory coordinates X in cell (one cell per channel), each row is a trajectory, missing values are nan
% dataset.Y: trajectory coordinates X in cell (one cell per channel), each row is a trajectory, missing values are nan
% dataset.Y: trajectory coordinates X in cell (one cell per channel), each row is a trajectory, missing values are nan
% dataset.interaction_data contains information about interactions 

% moviename : base name of the file that will be generated and contains 

% to_be_computed : structure that specify which quantities should be
% computed (see Examples, defined in the 'parameter' files).
% example:
% to_be_computed=struct;
% to_be_computed.exist=1;
% to_be_computed.inCCP=1;
% to_be_computed.inActin=0;
% to_be_computed.trapped=1;
% to_be_computed.inCCP_trapping_corrected=1;
% to_be_computed.inActin_trapping_corrected=0;
% to_be_computed.interaction=1;
% to_be_computed.binary_vec=1;
% to_be_computed.on_membrane=1;
% to_be_computed.state_number=1;
% to_be_computed.markov_chain_forward=1;
% to_be_computed.markov_chain_forward_supdur=1;

% list_channel : list of channel for each movie on which to compute the states list_channel=[1,2];

% list_interaction_couple : specifies interaction between which channels [1,2] for interaction between channel 1 and 2, can also be [1,2;1,3] for interaction between 1and 2 and between 1 and 3.
% must be [] if no interactions are calculated

% parameter : structure containing all parameters to be used in format parameter.***

% save_path : path to the folder where list of states will be stored 

%%

channel={'-C1', '-C2', '-C3', '-C4'};
list_of_file={'_gui2_steps.mat','_intmatrix_0pShiftX-0pShiftY-0fShiftT.mat'};
list_not_working={};

%%

disp(['Loading ',moviename])

%% remove previous error file for this movie
if isfile([save_path,filesep,moviename,'_ERROR_list_state.mat'])
    delete([save_path,filesep,moviename,'_ERROR_list_state.mat'])
end
error_movie={};%create empty structure for possible errors
list_state_ch=cell(1,max(list_channel));

for ch=list_channel
    % list_state_ch{ch}=struct;
    movie_name=[moviename,channel{1,ch}];
    %% check that file list_state exists, if it doesn't, the cript creates an empty structure and save it
    if ~isfile([save_path,filesep,movie_name,'_list_state.mat'])
        disp([movie_name, ' creating list of states'])
        list_state=struct;
        
    else
        load([save_path,filesep,movie_name,'_list_state.mat'])
        if exist('list_state','var')==1
            disp([movie_name, ' loaded list of states'])
        else
            disp([movie_name, ' creating list of states'])
            list_state=struct;
        end
    end
    list_state_ch{ch}=list_state;
    clearvars list_state
end
for ch=list_channel
    movie_name=[moviename,channel{1,ch}];
    %% look at the states to be computed
        
        if isfield(list_state_ch{ch},'exist')==0
            test.exist=1;
        else
            test.exist=0;
        end
        
        if isfield(list_state_ch{ch},'inCCP')==0
            test.inCCP=1;
        else
            test.inCCP=0;
        end
        
        if isfield(list_state_ch{ch},'inActin')==0
            test.inActin=1;
        else
            test.inActin=0;
        end
        
        if isfield(list_state_ch{ch},'trapped')==0
            test.trapped=1;
        else
            test.trapped=0;
        end
        %%
        if isfield(list_state_ch{ch},'inCCP_trapping_corrected')==0
            test.inCCP_trapping_corrected=1;
        else
            if list_state_ch{ch}.inCCP_trapping_corrected==1
                test.inCCP_trapping_corrected=0;
            else
                test.inCCP_trapping_corrected=1;
            end
        end
        %%
        if isfield(list_state_ch{ch},'inActin_trapping_corrected')==0
            test.inActin_trapping_corrected=1;
        else
            if list_state_ch{ch}.inActin_trapping_corrected==1
                test.inActin_trapping_corrected=0;
            else
                test.inActin_trapping_corrected=1;
            end
        end
        %%
        if isfield(list_state_ch{ch},'interaction')==0 && ~isempty(list_interaction_couple)
            test.interaction=1;
        else
            test.interaction=0;
        end
        
        if isfield(list_state_ch{ch},'on_membrane')==0
            test.on_membrane=1;
        else
            test.on_membrane=0;
        end
        if isfield(list_state_ch{ch},'binary_vec')==0
            test.binary_vec=1;
        else
            test.binary_vec=0;
        end
        
        clearvars list_state

    %% if the field 'exist' does not exist in the list of state, it is computed
    if to_be_computed.exist==1
        if test.exist==1
            try
                disp([movie_name, ' calculating list of exist'])
                [list_state_ch{ch}]=add_list_exist(dataset.X{ch},dataset.Y{ch},list_state_ch{ch},parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='exist';
                list_not_working{size(list_not_working,1),3}=e;
                
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='exist';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of exist'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    
    %% if the field 'inCCP' does not exist in the list of state, it is computed
    if to_be_computed.inCCP==1
        if test.inCCP==1
            try
                disp([movie_name, ' calculating list of inCCP'])
                [list_state_ch{ch}]=add_list_inmask(dataset.X{ch},dataset.Y{ch},dataset.stack_binary_CCP,'inCCP',list_state_ch{ch});
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='inCCP';
                list_not_working{size(list_not_working,1),3}=e;
                
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='inCCP';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of inCCP'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the field 'inActin' does not exist in the list of state, it is computed
    if to_be_computed.inActin==1
        if test.inActin==1
            try
                disp([movie_name, ' calculating list of inActin'])
                [list_state_ch{ch}]=add_list_inmask(dataset.X{ch},dataset.Y{ch},dataset.stack_binary_Actin,'inActin',list_state_ch{ch},parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='inActin';
                list_not_working{size(list_not_working,1),3}=e;
                
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='inActin';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of inActin'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the field trapped' does not exist in the list of state, it is computed
    if to_be_computed.trapped==1
        if test.trapped==1
            try
                disp([movie_name, ' calculating list of istrapped'])
                [list_state_ch{ch}]=add_list_istrapped(dataset.X{ch},dataset.Y{ch},list_state_ch{ch},parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='trapped';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='trapped';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of istrapped'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the field 'CCP_trapping_corrected' does not exist in the list of state, it is computed
    if to_be_computed.inCCP_trapping_corrected==1
        if test.inCCP_trapping_corrected==1
            try
                disp([movie_name, ' calculating CCP_trapping_correction'])
                [list_state_ch{ch}]=add_mask_trapping_corrected(list_state_ch{ch},'inCCP');
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='CCP_trapping_correction';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='CCP_trapping_correction';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of CCP_trapping_correction'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the field 'inActin_trapping_corrected' does not exist in the list of state, it is computed
    if to_be_computed.inActin_trapping_corrected==1
        if test.inActin_trapping_corrected==1
            try
                disp([movie_name, ' calculating inActin_trapping_correction'])
                [list_state_ch{ch}]=add_mask_trapping_corrected(list_state_ch{ch},'inActin');
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='inActin_trapping_correction';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='inActin_trapping_correction';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of inActin_trapping_correction'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the field 'interaction' does not exist in the list of state, it is computed
    if to_be_computed.interaction==1
        if test.interaction==1
            for n_couple=1:size(list_interaction_couple,1)
                try
                    disp([movie_name, ' calculating list of interactions'])
                    [list_state_ch{ch}]=add_list_interaction(dataset.X,dataset.interaction_data,ch,list_interaction_couple(n_couple,:),list_state_ch{ch});
                catch e
                    list_not_working{size(list_not_working,1)+1,1}=movie_name;
                    list_not_working{size(list_not_working,1),2}=['interaction C',num2str(list_interaction_couple(n_couple,1)),' C',num2str(list_interaction_couple(n_couple,2))];
                    list_not_working{size(list_not_working,1),3}=e;
                    error_movie{size(error_movie,1)+1,1}=movie_name;
                    error_movie{size(error_movie,1),2}='interaction';
                    error_movie{size(error_movie,1),3}=e;
                    disp([movie_name, ' error with calculation of interaction'])
                    disp(e)
                    disp([e.stack(1).name,' line ',num2str(e.stack(1).line)])
                    return
                end
            end
        end
    end
    
    %% if the field 'flipped' does not exist in the list of state, it is computed
    if to_be_computed.binary_vec==1
        if test.binary_vec==1
            try
                disp([movie_name, ' combining states'])
                [list_state_ch{ch}]=add_list_binary_vec(list_state_ch{ch},list_interaction_couple,parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='flipped';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='flipped';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of flipped'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %% if the file 'on_membrane' does not exist in the list of state, it is computed
    if to_be_computed.on_membrane==1
        if test.on_membrane==1
            try
                disp([movie_name, ' calculating list of on_membrane'])
                [list_state_ch{ch}]=add_list_onmembrane(list_state_ch{ch},list_interaction_couple,parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='on_membrane';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='on_membrane';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of on_membrane'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %  end
end

for ch=list_channel
    
    
%     if  test.forced==0
        
        
        if isfield(list_state_ch{ch},'state_number')==0
            test.state_number=1;
        else
            test.state_number=0;
        end
        if isfield(list_state_ch{ch},'markov_chain_forward')==0
            test.markov_chain_forward=1;
        else
            test.markov_chain_forward=0;
        end
        if isfield(list_state_ch{ch},'markov_chain_forward_supdur')==0
            test.markov_chain_forward_supdur=1;
        else
            test.markov_chain_forward_supdur=0;
        end
        %   clearvars list_state
        
%     end
    %% if the field 'state_number' does not exist in the list of state, it is computed
    if to_be_computed.state_number==1
        if test.state_number==1
            try
                disp([movie_name, ' calculating list of state_numbers'])
                [list_state_ch{ch}]=add_list_state_number(list_state_ch,ch,list_interaction_couple,parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='state_number';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='state_number';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of state_number'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
            
        end
    end
    
    %% add forward markov chain
    if to_be_computed.markov_chain_forward==1
        if test.markov_chain_forward==1
            try
                disp([movie_name, ' calculating markov_chain_forward'])
                [list_state_ch{ch}]=add_markov_chain_forward(list_state_ch{ch},parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='markov_chain_forward';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='markov_chain_forward';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of markov_chain_forward'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    if to_be_computed.markov_chain_forward_supdur==1
        
        if test.markov_chain_forward_supdur==1
            try
                disp([movie_name, ' calculating markov_chain_forward_supdur'])
                [list_state_ch{ch}]=add_markov_chain_forward_supdur(list_state_ch{ch},parameter);
            catch e
                list_not_working{size(list_not_working,1)+1,1}=movie_name;
                list_not_working{size(list_not_working,1),2}='markov_chain_forward_supdur';
                list_not_working{size(list_not_working,1),3}=e;
                error_movie{size(error_movie,1)+1,1}=movie_name;
                error_movie{size(error_movie,1),2}='markov_chain_forward_supdur';
                error_movie{size(error_movie,1),3}=e;
                disp([movie_name, ' error with calculation of markov_chain_forward_supdur'])
                disp(e)
                save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
                return
            end
        end
    end
    
    %  end
end
for ch=list_channel
    list_state=list_state_ch{ch};
    save([save_path,filesep,moviename,channel{1,ch},'_list_state.mat'],'list_state','-v7.3');
end
if ~isempty(error_movie)==1
    save([save_path,filesep,moviename,'_ERROR_list_state.mat'],'error_movie')
end
end

