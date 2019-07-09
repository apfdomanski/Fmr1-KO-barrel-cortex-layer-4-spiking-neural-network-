% Simulation and analysis code to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
%
% Reference: 
% Domanski,Booker,Wyllie,Isaac,Kind (2018)
% "Cellular and Synaptic Compensations Limit Circuit Disruption in Fmr1-KO Mouse but Fail to Prevent Deficits in Informaiton Processing"
%  BiorXiv preprint: https://doi.org/10.1101/403725
%
% Author: Aleksander PF Domanski 2015-2019 University of Bristol, UK aleks.domanski@bristol.ac.uk
%
% Copyright: (C) Aleksander PF Domanski 2019 University of Bristol, UK
%
% Usage:
% This script recapitulates the simulation results for Fig 9 in the above paper.
% - A recurrent spiking network model of barrel cortex Layer 4 is stimulated with a model thalamocortical bulk input.
% - Network architecture features sparse, randomly-connected Excitatory and Inhibiory leaky Integrate-and-Fire neurons
% - All synapses are modelled as sum-of-exponential currents with short-term plasticity and delayed transmission.
% - Synaptic, network connectivity and intrinsic parameters are parameterised by slice electrophysiology recordings detailed in Figs. 1-9 in the paper.
% - Parameter values are specified in input arguments and code for L4sim_DesignNetwork.m
% - This loop simulates ten repetitive trials of model TC stimulation
%   for ten random network seeds from both genotypes.  
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% http://www.gnu.org/copyleft/gpl.html

%% Startup, instantiate parallel pool
clear all; close all
if isempty(gcp('nocreate'))
    parpool local
end
if ispc
    filePath = ['C:\AD_data\ModellingResults' filesep datestr(datetime('now'), 'mm-dd-yyyy HH-MM-SS') ];
else ismac
    filePath = ['/Users/aleksanderdomanski/Documents/AD_DATA/ModellingResults/' datestr(datetime('now'), 'mm-dd-yyyy HH-MM-SS') ];
end
mkdir(filePath)
%% loop across genotypes
GenotypeList = {'WT','KO'};
no_models=10;
no_trials=10;
stimFreqs=[5 10 20 50];
%Simulation loop hierarchy:
% Genotype > Random model instantiation > Stimulation frequency > Trial no.
for iGenotype = 1:length(GenotypeList) 
    Genotype = GenotypeList{iGenotype};
    %% precalculate models
    switch Genotype
        case 'WT'
            for model_id = 1:no_models
                Marray{model_id}    = L4sim_DesignNetwork(800,150,...  % number of Ex, In neurons
                    0.7,2000,...               % TC->In STP parameters: dep. per stim and recovery time constant
                    0.7,500,...                % TC->Ex STP parameters: dep. per stim and recovery time constant
                    0.6,500,...                % In->Ex STP parameters: dep. per stim and recovery time constant
                    0.45,500,...               % Ex->In STP parameters: dep. per stim and recovery time constant
                    0.6,500,...                % In->In STP parameters: dep. per stim and recovery time constant
                    0.6,2000,...               % Ex->Ex STP parameters: dep. per stim and recovery time constant
                    [0.18,0.0001,0.0001],...     % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]
                    [0.5, 0.0001],...           % [p(Ex>In), mean strength]
                    [0.55,0.0002],...           % [p(In>Ex), mean strength]
                    [0.5, 0.0005],...           % [p(In>In), mean strength]
                    [0.005,0.025,0.005],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                    0,'WT');
            end
        case 'KO'
            for model_id = 1:no_models
                Marray{model_id} = Iz_design_con_matrix_Fmr1_COBN(800,150,...  % number of Ex, In neurons
                    0.2,8000,...                 % TC->In STP parameters: dep. per stim and recovery time constant
                    0.5,500,...                 % TC->Ex STP parameters: dep. per stim and recovery time constant
                    0.6,3000,...                 % In->Ex STP parameters: dep. per stim and recovery time constant
                    0.2,500,...                 % Ex->In STP parameters: dep. per stim and recovery time constant
                    0.6,3000,...                 % In->In STP parameters: dep. per stim and recovery time constant
                    0.6,5e3,...                 % Ex->Ex STP parameters: dep. per stim and recovery time constant
                    [0.18,0.0001,0.0001],...      % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]
                    [0.2,0.00015],...             % [p(Ex>In), mean strength]
                    [0.3,0.0002],...             % [p(In>Ex), mean strength]
                    [0.3,0.0005],...             % [p(In>In), mean strength]
                    [0.005,0.025,0.005],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                    0,'KO');
            end
    end
    
    parfor model_id = 1:no_models
        tic
        for iStim=1:length(stimFreqs)
            
            stimPattern = [5,  stimFreqs(iStim)];
            M    = L4sim_MakePulseInput(Marray{model_id},1000,0.5,[1 1],stimPattern,0);
            %%
            for trial_id=1:no_trials
                %% Runtime
                results=struct;
                
                %Produce trial-trial variability
                M.p.vr = M.p.vr + 0.5*randn(M.p.Ne+M.p.Ni,1); % randomly seed Vm to introduce trial-trial variability
                r=randperm(M.p.Ne)'; % Shuffle cells receiving TC stim
                M.net.S_ext(1:M.p.Ne)=M.net.S_ext(r);
                M.net.S_ext_NMDA(1:M.p.Ne)=M.net.S_ext_NMDA(r);
                
                fprintf('Running %s model no. %d, stim pattern %d x %dHz, trial no. %d...\n',...
                    Genotype,model_id,stimPattern(1),stimPattern(2),trial_id)
                
                results  = L4sim_RunModel(M,0,1);
                results  = L4sim_Analyse(M,results);
                % L4sim_PlotResults(model,results);
                
                Fname = sprintf('%s_%d_%d_%d.mat',Genotype,model_id,stimFreqs(iStim),trial_id)
                
                save_parfor([filePath filesep  Fname],'M')
                save_parfor([filePath filesep  Fname],'results')
                save_parfor([filePath filesep  Fname],'no_models')
                save_parfor([filePath filesep  Fname],'model_id')
                save_parfor([filePath filesep  Fname],'stimFreqs')
                save_parfor([filePath filesep  Fname],'iStim')
                save_parfor([filePath filesep  Fname],'no_trials')
                save_parfor([filePath filesep  Fname],'trial_id')
                
                %close all
            end
        end
    end
    
    clear trial_id model_id iStim Marray
end
delete(gcp)