% Simulated example short-term plasticity of model synapse for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% This scripts plots the short-term depression (and recovery from
% depression) at synapses between neurons in the model network
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
clear all; close all
%% WT model

M  = L4sim_DesignNetwork(10,10,...  % number of Ex, In neurons
                            0.6,5500,...               % TC->I  n STP parameters: dep. per stim and recovery time constant
                            0.67,5500,...                % TC->Ex STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->Ex STP parameters: dep. per stim and recovery time constant
                            0.6,3500,...               % Ex->In STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->In STP parameters: dep. per stim and recovery time constant
                            0.6,2000,...               % Ex->Ex STP parameters: dep. per stim and recovery time constant
                            [0.18,0.00013,0.00013],...     % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]           
                            [0.5, 0.00013],...           % [p(Ex>In), mean strength]           
                            [0.55,0.0002],...           % [p(In>Ex), mean strength]           
                            [0.5, 0.00020],...           % [p(In>In), mean strength]  
                            [0.00,0.25,0.00],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                            0,'WT');  


M        = L4sim_MakePulseInput(M,250,0.5,[1 1],[5 50],1);
results  = L4sim_RunModel(M,0,1); % (Model, online plot type, turn on depression Y/N)
%%
figure
subplot(2,1,1);hold on
plot(100*results.I_in_total(10,:)','b','LineWidth',1.5)
title('Synaptic current')
ylabel('Synaptic IPSC (nA)')
legend('WT Inhibitory synapse','Location','northeast');legend boxoff
subplot(2,1,2);hold on
plot(results.V_out(10,:)','b','LineWidth',1.5)
title('Synaptic potential')
ylabel('Synaptic IPSP (mV)')
xlabel('Time (ms)')
%% KO model
M = L4sim_DesignNetwork(10,10,...  % number of Ex, In neurons
                            0.4,5500,...                 % TC->In STP parameters: dep. per stim and recovery time constant
                            0.53,5500,...                 % TC->Ex STP parameters: dep. per stim and recovery time constant
                            0.6,3000,...                 % In->Ex STP parameters: dep. per stim and recovery time constant
                            0.2,500,...                 % Ex->In STP parameters: dep. per stim and recovery time constant
                            0.6,3000,...                 % In->In STP parameters: dep. per stim and recovery time constant
                            0.6,5e3,...                 % Ex->Ex STP parameters: dep. per stim and recovery time constant
                            [0.18,0.00012,0.00012],...      % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]           
                            [0.2,0.00015],...             % [p(Ex>In), mean strength]           
                            [0.3,0.0002],...             % [p(In>Ex), mean strength]           
                            [0.3,0.0005],...             % [p(In>In), mean strength]  
                            [0.00,0.1,0.0],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                            0,'KO'); 
                        
M        = L4sim_MakePulseInput(M,250,0.5,[1 1],[5 50],0);
results  = L4sim_RunModel(M,0,1); % (Model, online plot type, turn on depression Y/N)
%%
figure
subplot(2,1,1);hold on
plot(100*results.I_in_total(10,:)','r','LineWidth',1.5)
title('Synaptic current')
ylabel('Synaptic IPSC (nA)')
legend('KO Inhibitory synapse','Location','northeast');legend boxoff
subplot(2,1,2);hold on
plot(results.V_out(10,:)','r','LineWidth',1.5)
title('Synaptic potential')
ylabel('Synaptic IPSP (mV)')
xlabel('Time (ms)')         
                        