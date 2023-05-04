% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%% parameters
parameter=struct;
parameter.pixel_size=106.7;%nanometer
parameter.FrameNumber=1000;% number of frame to be analysed
parameter.Frameduration=0.03;% duration of a frame in second

parameter.min_traj_length=3;% minimum trajectory length to compute states

%%
% inMask
%                   'state_name','file suffix.format','position of the transformation piecewise linear param'(nan for no alignement)
parameter.mask_pairs={'inCCP','_stack_binary_mask_CCP.mat',3;...
'inActin','-actin_msk2.tif',nan};


%% interaction
parameter.list_interaction_couple=[1,2];
parameter.frame_number_interaction_analysis=1000;
parameter.searchRadius=150;%nm
parameter.searchRadius_pixel=parameter.searchRadius/parameter.pixel_size;%nm


%% parameter trapping detection
parameter.p_val_traj_type='fBm'; % can be 'Bm' or 'fBm'
parameter.T_mean=2; % must stay 2
parameter.sig_noise=0; % must stay 0
parameter.diag_percentile=10;% can be 0,5,10 or 50
parameter.nu=0.75; % can be 0.1,0.3,0.5,0.75,0.9,1; 0.75 is recommended
parameter.p_value=0.05; % can be any percentile (minimum 0.01)
parameter.list_mu=[1:0.5:2]; % maximum range ([0.5,1,1.5,2,2.5,3])

%% Complete list of states
% flipped :> CCP | istrapped | interaction
parameter.state_vec_fields={'inCCP','trapped','interaction','on_membrane'}; % which fields of list_states are used in which order for state number attribution
% dummy state for frame before molecule appears and after it disappears
parameter.state{1,1}="Absent";
parameter.state{1,2}={[-1;-1;-1;0]};
parameter.state{1,3}={[]};

% not trapped
parameter.state{2,1}="Free"; % if molecule is not trapped and not colocalising, or colocalising with a trapped molecule
parameter.state{2,2}={[0;0;0;1],[1;0;0;1],[0;0;1;1],[0;0;1;1],[1;0;1;1],[1;0;1;1]};
parameter.state{2,3}={[],[],[0;1;1;1],[1;1;1;1],[0;1;1;1],[1;1;1;1]};

parameter.state{3,1}="Co-diffusion"; % if free molecule colocalise with a free partner
parameter.state{3,2}={[0;0;1;1],[1;0;1;1],[1;0;1;1],[0;0;1;1]};
parameter.state{3,3}={[0;0;1;1],[0;0;1;1],[1;0;1;1],[1;0;1;1]};

% trapped outside CCP
parameter.state{4,1}="Confined";% if molecule is confined outside CCP and (is alone or colocalise with a free partner)
parameter.state{4,2}={[0;1;0;1],[0;1;1;1],[0;1;1;1]};
parameter.state{4,3}={[],[0;0;1;1],[1;0;1;1]};

parameter.state{5,1}="Co-confined"; % if molecule and its partner are both confined outside CCP
parameter.state{5,2}={[0;1;1;1],[0;1;1;1]};
parameter.state{5,3}={[0;1;1;1],[1;1;1;1]};

% trapped inside CCP
parameter.state{6,1}="Trapped alone in CCP";% if molecule is confined inside CCP and (is alone or colocalise with a free partner)
parameter.state{6,2}={[1;1;0;1],[1;1;1;1],[1;1;1;1]};
parameter.state{6,3}={[],[0;0;1;1],[1;0;1;1]};

parameter.state{7,1}="Co-trapped in CCP"; % if molecule and its partner are both trapped inside CCP
parameter.state{7,2}={[1;1;1;1],[1;1;1;1]};
parameter.state{7,3}={[0;1;1;1],[1;1;1;1]};

%% Inputs for state analysis

% which part of analysis to be computed
to_be_computed=struct;
to_be_computed.exist=1;
to_be_computed.inCCP=1;
to_be_computed.inActin=0;
to_be_computed.trapped=1;
to_be_computed.inCCP_trapping_corrected=1;
to_be_computed.inActin_trapping_corrected=0;
to_be_computed.interaction=1;
to_be_computed.binary_vec=1;
to_be_computed.on_membrane=1;
to_be_computed.state_number=1;
to_be_computed.markov_chain_forward=1;
to_be_computed.markov_chain_forward_supdur=1;

%% parameter for Markov Chain
parameter.markov_chain_Nmax=200;
parameter.sup_dur=10;% minimum duration of a state in markov_chain_supdur
