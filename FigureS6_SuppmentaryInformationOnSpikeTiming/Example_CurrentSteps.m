% Simulated example intrinsic excitability of model neurons for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% This scripts simulates the whole-cell current-injection results reported in the
% accompanying paper by looping across current step amplitudes and
% providing step current input.
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
deltaCurrent = 10;
steps = deltaCurrent:deltaCurrent:200;
stimAmp = steps*1e-3;
%% Run model
% Run WT model
output.noSpikes=[];
M = L4sim_DesignNetwork(100,100,...  % number of Ex, In neurons
                            0.6,5500,...               % TC->I  n STP parameters: dep. per stim and recovery time constant
                            0.67,5500,...                % TC->Ex STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->Ex STP parameters: dep. per stim and recovery time constant
                            0.6,3500,...               % Ex->In STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->In STP parameters: dep. per stim and recovery time constant
                            0.6,2000,...               % Ex->Ex STP parameters: dep. per stim and recovery time constant
                            [1e-10,0,0],...     % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]           
                            [1e-10,0],...           % [p(Ex>In), mean strength]           
                            [1e-10,0],...           % [p(In>Ex), mean strength]           
                            [1e-10, 0],...           % [p(In>In), mean strength]  
                            [0.006,0.025,0.0],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                            0,'WT');  

for iStim=1:length(stimAmp)
    M.p.vr=-70*ones(size(M.p.vr));
    results=struct;
    %%% Runtime
    M        = L4sim_MakeCurentInjectionInput(M,1000,1,[stimAmp(iStim) stimAmp(iStim)],500,0);
    results  = L4sim_RunCurrentInjection(M,0,0); % (Model, online plot type, turn on depression Y/N)
    output.firings{iStim}=results.firings;
    output.Vout{iStim}=results.V_out;
end
for iStim=1:length(stimAmp)
    if ~isempty(output.firings{iStim})  
        for iNeuron=1:M.p.Ne+M.p.Ni
            output.noSpikes(iNeuron,iStim) = sum(output.firings{iStim}(:,2)==iNeuron);
        end
    else
        output.noSpikes(1:M.p.Ne+M.p.Ni,iStim) = 0;
    end
end
WT=output;

% run KO model
output.noSpikes=[];
M = L4sim_DesignNetwork(1,1,...  % number of Ex, In neurons
                            0.6,5500,...               % TC->I  n STP parameters: dep. per stim and recovery time constant
                            0.67,5500,...                % TC->Ex STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->Ex STP parameters: dep. per stim and recovery time constant
                            0.6,3500,...               % Ex->In STP parameters: dep. per stim and recovery time constant
                            0.8,3000,...                % In->In STP parameters: dep. per stim and recovery time constant
                            0.6,2000,...               % Ex->Ex STP parameters: dep. per stim and recovery time constant
                            [1,0,0],...     % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]           
                            [1,0],...           % [p(Ex>In), mean strength]           
                            [1,0],...           % [p(In>Ex), mean strength]           
                            [1e-10, 0],...           % [p(In>In), mean strength]  
                            [0.006,0.025,0.0],...      % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                            0,'KO');  

for iStim=1:length(stimAmp)
    M.p.vr=-70*ones(size(M.p.vr));
    results=struct;
    %%% Runtime
    M        = L4sim_MakeCurentInjectionInput(M,1000,1,[stimAmp(iStim) stimAmp(iStim)],500,0);
    results  = L4sim_RunCurrentInjection(M,0,0); % (Model, online plot type, turn on depression Y/N)
    output.firings{iStim}=results.firings;
    output.Vout{iStim}=results.V_out;
end
for iStim=1:length(stimAmp)
    if ~isempty(output.firings{iStim})  
        for iNeuron=1:M.p.Ne+M.p.Ni
                output.noSpikes(iNeuron,iStim) = sum(output.firings{iStim}(:,2)==iNeuron);
        end
    else
        output.noSpikes(1:M.p.Ne+M.p.Ni,iStim) = 0;
    end
end
KO=output;
%% Plot F-I plot (no. spikes vs current injection)
WT.noSpikes_Exmean=mean(WT.noSpikes(1:M.p.Ne,:),1);
KO.noSpikes_Exmean=mean(KO.noSpikes(1:M.p.Ne,:),1);

WT.noSpikes_ExSEM=nansem(WT.noSpikes(1:M.p.Ne,:),1);
KO.noSpikes_ExSEM=nansem(KO.noSpikes(1:M.p.Ne,:),1);

WT.noSpikes_Inmean=mean(WT.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
KO.noSpikes_Inmean=mean(KO.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);

WT.noSpikes_InSEM=nansem(WT.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
KO.noSpikes_InSEM=nansem(KO.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);



figure; 
subplot(1,2,1);hold on

plot(steps,KO.noSpikes(1:M.p.Ne,:),'color',[0.9 0.6 0.6])
plot(steps,WT.noSpikes(1:M.p.Ne,:),'color',[0.6 0.6 0.9])

errorbar(steps,WT.noSpikes_Exmean,WT.noSpikes_ExSEM,'LineWidth',2,'color','b')
errorbar(steps,KO.noSpikes_Exmean,KO.noSpikes_ExSEM,'LineWidth',2,'color','r')
axis([0 max(stimAmp)*1e3 0 100])
xlabel('Current step (pA)')
ylabel('No. Spikes Fired /500ms')
title('Ex neurons')

subplot(1,2,2);hold on

plot(steps,KO.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),'color',[0.9 0.6 0.6])
plot(steps,WT.noSpikes(M.p.Ne+1:M.p.Ne+M.p.Ni,:),'color',[0.6 0.6 0.9])

errorbar(steps,WT.noSpikes_Inmean,WT.noSpikes_InSEM,'LineWidth',2,'color','b')
errorbar(steps,KO.noSpikes_Inmean,KO.noSpikes_InSEM,'LineWidth',2,'color','r')
axis([0 max(stimAmp)*1e3 0 100])
xlabel('Current step (pA)')
ylabel('No. Spikes Fired /500ms')
title('FS neurons')
legend('WT (Fmr1 +/Y) model','KO (Fmr1 -/Y model)','Location','northwest'); legend boxoff
%% Plot example Ex neuron

WTtemp=[];KOtemp=[];
for iStim=1:length(stimAmp)
    
    WTtemp=[WTtemp;WT.Vout{iStim}(1,:)+120*iStim ];
    KOtemp=[KOtemp;KO.Vout{iStim}(1,:)+120*iStim];
end
x=1:1000;
figure; hold on
plot(x,WTtemp','b','LineWidth',1)
plot(x+1000,KOtemp','r','LineWidth',1)
axis([0 2000 -70 Inf])
axis off
%% Plot example FS neuron

WTtemp=[];KOtemp=[];
for iStim=1:length(stimAmp)
    
    WTtemp=[WTtemp;WT.Vout{iStim}(M.p.Ne+1,:)+120*iStim ];
    KOtemp=[KOtemp;KO.Vout{iStim}(M.p.Ne+1,:)+120*iStim];
end
x=1:1000;
figure; hold on
plot(x,WTtemp','b','LineWidth',1)
plot(x+1000,KOtemp','r','LineWidth',1)
axis([0 2000 -70 Inf])
axis off