% Analysis code to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% This script recapitulates the group analysis of simulation results for 
% Fig 9B-D in the above paper and plots example firing rates, simulated
% membrane potential and dependence of input frequency on network firing
% rates.
%
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

if ispc
    pat= 'C:\AD_data\ModellingResults';
else ismac
    pat = '/Users/aleksanderdomanski/Documents/AD_DATA/ModellingResults/';
end

stimFreqs=[5 10 20 50];
noModels=5;
noStimConds=length(stimFreqs);
noTrials=10;
recFolderWT = uigetdir(pat,'Select WT folder');
recFolderKO = uigetdir(pat,'Select KO folder');

%% Load WT waveforms
cd(recFolderWT)
Highlight = 503;

for iModel=4%1:noModels%noModels;
    for iFreq = 1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('WT model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            resmatrix.WT(iModel).results{iFreq,iTrial}=results.V_out(Highlight,:);
        end
    end
end
%% Load KO waveforms
cd(recFolderKO)
Highlight = 417;%430;%

for iModel=4%1:noModels
    for iFreq = 1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('KO model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            resmatrix.KO(iModel).results{iFreq,iTrial}=results.V_out(Highlight,:);
        end
    end
end
%% plot WT waveforms
Highlight = 503;
nTotal=M.p.Ne;M.p.Ni;
for iModel= 4%1:noModels
    for iTrial = 1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2-1); hold on
            plot(resmatrix.WT(iModel).results{iFreq,iTrial}','b')
            axis([0 Inf -80 50]); axis off
        end
    end
end
%% plot KO waveforms
Highlight = 417;%430;%

nTotal=M.p.Ne;M.p.Ni;
for iModel= 4%1:noModels
    for iTrial = 1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2); hold on
            plot(resmatrix.KO(iModel).results{iFreq,iTrial}','r')
            axis([0 Inf -80 50])
        end
    end
end
%% plot both waveforms

for iModel= 4%1:noModels
    for iTrial = 2%1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2-1); hold on
            plot(resmatrix.WT(iModel).results{iFreq,iTrial}(1,:)','b')
            axis([0 Inf -80 50]); axis off
        end
    end
end

nTotal=M.p.Ne;M.p.Ni;
for iModel= 4%1:noModels
    for iTrial = 2%1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2); hold on
            plot(resmatrix.KO(iModel).results{iFreq,iTrial}(1,:)','r')
            axis([0 Inf -80 50]); axis off
        end
    end
end
%% load WT raster
recFolderWT = uigetdir(pat,'Select WT folder');
cd(recFolderWT)
for iModel=4%1:noModels%noModels;
    for iFreq = 1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('WT model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            resmatrix.WT(iModel).firings{iFreq,iTrial}=results.firings;
            resmatrix.WT(iModel).raster{iFreq,iTrial}=sparse(results.spike_raster);
        end
    end
end
%% plot WT raster

Ne=M.p.Ne;  Ni=M.p.Ni;  N_total= Ne+Ni; dt=0.5;
for iModel= 4%1:noModels
    for iTrial = 1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2-1); hold on
            temp=resmatrix.WT(iModel).firings{iFreq,iTrial};
            temp_E=temp(temp(:,2)<=Ne,:);
            temp_I=temp(temp(:,2)>Ne,:);
            plot(temp_E(:,1),temp_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',8);
            plot(temp_I(:,1),temp_I(:,2),'.','color',[0.6 0.6 0.9],'MarkerSize',8);
            area((0:dt:1000-dt),...
                10*sum(resmatrix.WT(iModel).raster{iFreq,iTrial}(1:Ne,:)),'FaceColor',[0.6 0.6 0.6])
            %             area((0:dt:1000-dt),...
            %                 10*sum(resmatrix.WT(iModel).raster{iFreq,iTrial}(Ne+1:N_total,:)),'FaceColor',[0.6 0.6 0.9])
            axis([0 1000 0 N_total])
            %             axis off
        end
    end
end
%% load KO raster
recFolderKO = uigetdir(pat,'Select KO folder');
cd(recFolderKO)
for iModel=4%1:noModels%noModels;
    for iFreq = 1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('KO model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            resmatrix.KO(iModel).firings{iFreq,iTrial}=results.firings;
            resmatrix.KO(iModel).raster{iFreq,iTrial}=sparse(results.spike_raster);
        end
    end
end
%% plot KO raster
Ne=M.p.Ne;  Ni=M.p.Ni;  N_total= Ne+Ni;
for iModel= 4%:noModels
    for iTrial = 1:noTrials
        for iFreq = 1:noStimConds
            subplot(noStimConds,2,(iFreq)*2); hold on
            temp=resmatrix.KO(iModel).firings{iFreq,iTrial};
            temp_E=temp(temp(:,2)<=Ne,:);
            temp_I=temp(temp(:,2)>Ne,:);
            plot(temp_E(:,1),temp_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',8);
            plot(temp_I(:,1),temp_I(:,2),'.','color',[0.6 0.6 0.9],'MarkerSize',8);
            area((0:dt:1000-dt),...
                10*(sum(resmatrix.KO(iModel).raster{iFreq,iTrial}(1:Ne,:))),'FaceColor',[0.6 0.6 0.6])
            area((0:dt:1000-dt),...
                10*sum(resmatrix.KO(iModel).raster{iFreq,iTrial}(Ne+1:N_total,:)),'FaceColor',[0.6 0.6 0.9],'EdgeColor',[0.6 0.6 0.9])
            axis([0 1000 0 N_total])
            %             axis off
        end
    end
end
%% Load specified WT waveforms
iModel = 4 ; iTrial=7;
cd(recFolderWT)
for iFreq = 1:noStimConds
    
    disp(sprintf('WT model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
    A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
    load(A(1).name);
    resmatrix.WT(iModel).results{iFreq,iTrial}=results.V_out;
end
%% Load specified KO waveforms
iModel = 4 ; iTrial=8;
cd(recFolderKO)
for iFreq = 1:noStimConds
    
    disp(sprintf('KO model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
    A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
    load(A(1).name);
    resmatrix.KO(iModel).results{iFreq,iTrial}=results.V_out;
end
%% Plot both waveforms - all cells
figure('color','w');

Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
Highlight = 503;
iModel = 4 ; iTrial=7;
for iFreq = 1:noStimConds
    subplot(noStimConds,2,(iFreq)*2-1); hold on
    plot((0:dt:1000-dt), resmatrix.WT(iModel).results{iFreq,iTrial}(1:Ne-1,:)','color',[0.6 0.6 0.6])
    plot((0:dt:1000-dt),resmatrix.WT(iModel).results{iFreq,iTrial}(Highlight,:)','color',[0.2 0.2 0.9],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
Highlight = 417;%430;%
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
iModel = 4 ; iTrial=8;

for iFreq = 1:noStimConds
    subplot(noStimConds,2,(iFreq)*2); hold on
    plot((0:dt:1000-dt),resmatrix.KO(iModel).results{iFreq,iTrial}(1:Ne-1,:)','color',[0.6 0.6 0.6])
    plot((0:dt:1000-dt),resmatrix.KO(iModel).results{iFreq,iTrial}(Highlight,:)','color',[0.9 0.2 0.2],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
plot([750 850],[-20 -20],'k','LineWidth',1.5)
plot([750 750],[-20 0],'k','LineWidth',1.5)
text(680, -80,'100ms/20mV')
%% Plot both waveforms - one cell
figure('color','w');

Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
Highlight = 503;
iModel = 4 ; iTrial=7;
for iFreq = 1:noStimConds
    %     subplot(noStimConds,2,(iFreq)*2-1); hold on
    subplot(noStimConds,1,iFreq); hold on
    plot((0:dt:1000-dt),resmatrix.WT(iModel).results{iFreq,iTrial}(Highlight,:)','color',[0.5 0.5 0.9],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
Highlight = 417;%430;%
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
iModel = 4 ; iTrial=8;

for iFreq = 1:noStimConds
    %     subplot(noStimConds,2,(iFreq)*2); hold on
    subplot(noStimConds,1,iFreq); hold on
    plot((0:dt:1000-dt),resmatrix.KO(iModel).results{iFreq,iTrial}(Highlight,:)','color',[0.9 0.5 0.5],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
plot([750 850],[-20 -20],'k','LineWidth',1.5)
plot([750 750],[-20 30],'k','LineWidth',1.5)
text(680, -80,'100ms/50mV')
%% Plot both waveforms - one cell offset
figure('color','w');

Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
Highlight = 503;
iModel = 4 ; iTrial=7;
for iFreq = 1:noStimConds
    %     subplot(noStimConds,2,(iFreq)*2-1); hold on
    subplot(noStimConds,1,iFreq); hold on
    temp=resmatrix.WT(iModel).results{iFreq,iTrial}(Highlight,:)'; %temp(temp> -20)=-20;
    plot((0:dt:1000-dt),temp,'color',[0.2 0.2 0.9],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
Highlight = 417;%430;%
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
iModel = 4 ; iTrial=8;

for iFreq = 1:noStimConds
    %     subplot(noStimConds,2,(iFreq)*2); hold on
    subplot(noStimConds,1,iFreq); hold on
    temp=resmatrix.KO(iModel).results{iFreq,iTrial}(Highlight,:)'; %temp(temp> -20)=-20;
    plot((0:dt:1000-dt),temp,'color',[0.9 0.2 0.2],'LineWidth',1.5)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-90 -90 -90 -90 -90],'.k')
    axis([0 1000 -90 100])
    axis off
end
plot([750 850],[10 10],'k','LineWidth',1.5)
plot([750 750],[10 60],'k','LineWidth',1.5)
text(680, -80,'100ms/50mV')
%% Plot both waveforms - 5Hz only
figure('color','w');
hold on
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;


Highlight = 417;%430;%
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
iModel = 4 ; iTrial=8;

for iFreq = 1%:noStimConds
    %     subplot(noStimConds,2,(iFreq)*2); hold on
    %     subplot(noStimConds,1,iFreq); hold on
    temp=resmatrix.KO(iModel).results{iFreq,iTrial}(Highlight,:)'; temp(temp> -20)=-30;
    temp=temp-5
    plot((0:dt:1000-dt),temp,'color',[0.9 0.2 0.2],'LineWidth',2)
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-75 -75 -75 -75 -75],'.k')
    axis([0 1000 -90 0])
    axis off
end
Highlight = 503;
iModel = 4 ; iTrial=7; iFreq = 1%:noStimConds
temp=resmatrix.WT(iModel).results{iFreq,iTrial}(Highlight,:)'; %temp(temp> -20)=-20;
plot((0:dt:1000-dt),temp,'color',[0.2 0.2 0.9],'LineWidth',2)


plot([750 850],[-40 -40],'k','LineWidth',1.5)
plot([750 750],[-40 -30],'k','LineWidth',1.5)
text(680, -80,'100ms/10mV')

%% plot both rasters
figure('color','w');
markerSize=4;
Ne=M.p.Ne;  Ni=M.p.Ni;
Ne=800; Ni = 200; N_total= Ne+Ni; dt=0.5;
iModel = 4 ; iTrial=7;

for iFreq = 1:noStimConds
    subplot(noStimConds,2,(iFreq)*2-1); hold on
    temp=resmatrix.WT(iModel).firings{iFreq,iTrial};
    temp_E=temp(temp(:,2)<=Ne,:);
    temp_I=temp(temp(:,2)>Ne,:);
    plot(temp_E(:,1),temp_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',markerSize);
    plot(temp_I(:,1),temp_I(:,2),'.','color',[0.6 0.6 0.9],'MarkerSize',markerSize);
    area((0:dt:1000-dt),...
        10*(sum(resmatrix.WT(iModel).raster{iFreq,iTrial}(1:Ne,:))),'FaceColor',[0.6 0.6 0.6])
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-50 -50 -50 -50 -50],'.k')
    if iFreq == 1
        axis([00 1000 -50 N_total])
    else
        axis([00 500 -50 N_total])
    end
    axis off
    
end
iModel = 4 ; iTrial=8;
for iFreq = 1:noStimConds
    subplot(noStimConds,2,(iFreq)*2); hold on
    temp=resmatrix.KO(iModel).firings{iFreq,iTrial};
    temp_E=temp(temp(:,2)<=Ne,:);
    temp_I=temp(temp(:,2)>Ne,:);
    plot(temp_E(:,1),temp_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',markerSize);
    plot(temp_I(:,1),temp_I(:,2),'.','color',[0.6 0.6 0.9],'MarkerSize',markerSize);
    
    temp_E=10*(sum(resmatrix.KO(iModel).raster{iFreq,iTrial}(1:Ne,:)));
    temp_I=10*sum(resmatrix.KO(iModel).raster{iFreq,iTrial}(Ne+1:N_total,:));
    area((0:dt:1000-dt),temp_E,'FaceColor',[0.6 0.6 0.6])
    ISI=1000*(stimFreqs(iFreq)^-1);
    scatter([10,...
        10+ISI,...
        10+2*ISI,...
        10+3*ISI,...
        10+4*ISI], [-50 -50 -50 -50 -50],'.k')
    if iFreq == 1
        axis([00 1000 -50 N_total])
        plot([900 1000],[100 100],'k','LineWidth',1.5)
        text(800, 230,'100ms')
        
        
    else
        axis([00 500 -50 N_total])
    end
    axis off
end
plot([350 450],[100 100],'k','LineWidth',1.5)
plot([350 350],[100 500],'k','LineWidth',1.5)
text(380, 230,'50% co-active / 100ms')
% plot([350 350],[100 200],'k','LineWidth',1.5)
%% Total spike output /trial
analysis =struct;
for iModel=1:noModels%noModels;
    for iFreq = 1:noStimConds
        for iTrial = 1:noTrials
            temp=full(resmatrix.WT(iModel).raster{iFreq,iTrial});
            analysis.WT.sumSpikes{iFreq}(iModel,iTrial)= sum(sum(temp(1:Ne,:)'));
            analysis.WT.meanSpikes{iFreq}(iModel,iTrial)= mean(sum(temp(1:Ne,:)'));
            
            temp=full(resmatrix.KO(iModel).raster{iFreq,iTrial});
            analysis.KO.sumSpikes{iFreq}(iModel,iTrial)= sum(sum(temp(1:Ne,:)'));
            analysis.KO.meanSpikes{iFreq}(iModel,iTrial)= mean(sum(temp(1:Ne,:)'));
        end
    end
end

analysis.WT.sumSpikesMean=cellfun(@mean,analysis.WT.sumSpikes,'UniformOutput',false)
analysis.WT.sumSpikesMean=cell2mat(cellfun(@mean,analysis.WT.sumSpikesMean,'UniformOutput',false))

analysis.WT.sumSpikesSD=cellfun(@nansem,analysis.WT.sumSpikes,'UniformOutput',false)
analysis.WT.sumSpikesSD=cell2mat(cellfun(@mean,analysis.WT.sumSpikesSD,'UniformOutput',false))

analysis.KO.sumSpikesMean=cellfun(@mean,analysis.KO.sumSpikes,'UniformOutput',false)
analysis.KO.sumSpikesMean=cell2mat(cellfun(@mean,analysis.KO.sumSpikesMean,'UniformOutput',false))

analysis.KO.sumSpikesSD=cellfun(@nansem,analysis.KO.sumSpikes,'UniformOutput',false)
analysis.KO.sumSpikesSD=cell2mat(cellfun(@mean,analysis.KO.sumSpikesSD,'UniformOutput',false))

figure; hold on
plot([1 2 3 4],analysis.WT.sumSpikesMean,'-ob','LineWidth',1.5,'MarkerFaceColor','b')
errorbar([1 2 3 4],analysis.WT.sumSpikesMean,analysis.WT.sumSpikesSD,'b','LineWidth',1.5  )

plot([1 2 3 4],analysis.KO.sumSpikesMean,'-or','LineWidth',1.5,'MarkerFaceColor','r')
errorbar([1 2 3 4],analysis.KO.sumSpikesMean,analysis.KO.sumSpikesSD,'r','LineWidth',1.5  )
set(gca,'YScale','log')
set(gca,'Xtick',[1 2 3 4],...
    'XTickLabel',{'5' '10' '20' '50'})
%     title('Network output')
xlabel('Input frequency (Hz)')
ylabel('Network spike count /trial')
%% Spike count Histogram
analysis.WT.SpikeHist=[];
analysis.KO.SpikeHist=[];
for iFreq = 1:noStimConds
    for iModel=1:noModels%noModels;
        for iTrial = 1:noTrials
            temp=full(resmatrix.WT(iModel).raster{iFreq,iTrial});
            analysis.WT.SpikeHist{iFreq,iModel}(iTrial,:)= sum(temp(1:Ne,:));
            
            temp=full(resmatrix.KO(iModel).raster{iFreq,iTrial});
            analysis.KO.SpikeHist{iFreq,iModel}(iTrial,:)= sum(temp(1:Ne,:));
        end
    end
    temp = analysis.WT.SpikeHist(iFreq,:)';
    analysis.WT.SpikeHistMean{iFreq}=cell2mat(cellfun(@nanmax,temp,'UniformOutput',false));
    analysis.WT.SpikeHistSEM{iFreq}=cell2mat(cellfun(@nansem,temp,'UniformOutput',false));
    
    
    temp = analysis.KO.SpikeHist(iFreq,:)';
    analysis.KO.SpikeHistMean{iFreq}=cell2mat(cellfun(@nanmax,temp,'UniformOutput',false));
    analysis.KO.SpikeHistSEM{iFreq}=cell2mat(cellfun(@nansem,temp,'UniformOutput',false));
end

analysis.WT.SpikeHistGrandMean = cellfun(@nanmean,analysis.WT.SpikeHistMean,'UniformOutput',false);
analysis.WT.SpikeHistGrandSEM = cellfun(@nansem,analysis.WT.SpikeHistMean,'UniformOutput',false);
analysis.KO.SpikeHistGrandMean = cellfun(@nanmean,analysis.KO.SpikeHistMean,'UniformOutput',false);
analysis.KO.SpikeHistGrandSEM = cellfun(@nansem,analysis.KO.SpikeHistMean,'UniformOutput',false);
%% grand mean
figure;
hold on
iFreq=4;
area((0:dt:1000-dt),100*analysis.KO.SpikeHistGrandMean{iFreq}/Ne,'FaceColor',[0.9 0.6 0.6],'EdgeColor',[1 0 0],'FaceAlpha',1)
area((0:dt:1000-dt),100*analysis.WT.SpikeHistGrandMean{iFreq}/Ne,'FaceColor',[0.6 0.6 0.9],'EdgeColor',[0 0 1],'FaceAlpha',0.8)
ISI=1000*(stimFreqs(iFreq)^-1);
scatter([10,...
    10+ISI,...
    10+2*ISI,...
    10+3*ISI,...
    10+4*ISI], [-1 -1 -1 -1 -1],'.k')
axis([0 500 -2 20])
plot([350 450],[5 5],'k','LineWidth',1.5)
plot([350 350],[5 10],'k','LineWidth',1.5)
text(360, 7,{'10% Synchrony'; '100ms'})
axis off
%% example repeats
figure;
hold on
iFreq=4; iModel=4;
temp = 100*analysis.WT.SpikeHist{iFreq,iModel}/Ne+18;
plot((0:dt:1000-dt),temp(:,1:end),'color',[0.3 0.3 0.9],'LineWidth',1.5)

iFreq=4; iModel=4;
temp = 100*analysis.KO.SpikeHist{iFreq,iModel}/Ne;
plot((0:dt:1000-dt),temp(:,1:end),'color',[0.9 0.3 0.3],'LineWidth',1.5)
ISI=1000*(stimFreqs(iFreq)^-1);
scatter([10,...
    10+ISI,...
    10+2*ISI,...
    10+3*ISI,...
    10+4*ISI], [-1 -1 -1 -1 -1]+18,300,'.k')
scatter([10,...
    10+ISI,...
    10+2*ISI,...
    10+3*ISI,...
    10+4*ISI], [-1 -1 -1 -1 -1],300,'.k')
axis([50 150 -2 Inf])
axis off

% plot([10 20],[5 5],'k','LineWidth',1.5)
plot([130 130],[5 10],'k','LineWidth',2)
text(135, 7,{'5% Synchrony'})
%% spike stats
iFreq=4;
ISIbins=0:1:200;
analysis.WT.FiringSpan = [];
analysis.WT.FirstSpike = [];
analysis.WT.NoSpikes   = [];
analysis.WT.ISIhist    = [];
for iModel = 1:noModels
    for iTrial = 1:noTrials
        temp=resmatrix.WT(iModel).firings{iFreq,iTrial};
        for iNeuron= 1:Ne
            if length(temp(temp(:,2)==iNeuron,1))>1 %lots of spikes
                analysis.WT.FiringSpan{iModel,iTrial}(iNeuron,:) = [min(temp(temp(:,2)==iNeuron,1)), max(temp(temp(:,2)==iNeuron,1))];
                analysis.WT.FirstSpike{iModel}(iTrial,iNeuron)   = min(temp(temp(:,2)==iNeuron,1));
                analysis.WT.NoSpikes{iModel}(iTrial,iNeuron)     = length(temp(temp(:,2)==iNeuron,1));
                analysis.WT.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = histc(diff(temp(temp(:,2)==iNeuron,1)),ISIbins);
                
            elseif length(temp(temp(:,2)==iNeuron,1))==1 % 1 spike
                analysis.WT.FiringSpan{iModel,iTrial}(iNeuron,:) = [min(temp(temp(:,2)==iNeuron,1)), NaN];
                analysis.WT.FirstSpike{iModel}(iTrial,iNeuron)   = min(temp(temp(:,2)==iNeuron,1));
                analysis.WT.NoSpikes{iModel}(iTrial,iNeuron)     = 1;
                analysis.WT.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = nan(size(ISIbins));
                
                
            else % no spikes
                analysis.WT.FiringSpan{iModel,iTrial}(iNeuron,:) = [NaN,NaN];
                analysis.WT.FirstSpike{iModel}(iTrial,iNeuron)   = NaN;
                analysis.WT.NoSpikes{iModel}(iTrial,iNeuron)     = 0;
                analysis.WT.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = nan(size(ISIbins));
            end
        end
    end
end
analysis.WT.FirstSpike_mean = cellfun(@nanmean,analysis.WT.FirstSpike,'UniformOutput',false);
analysis.WT.FirstSpike_mean = cell2mat(cellfun(@transpose,analysis.WT.FirstSpike_mean,'UniformOutput',false));

analysis.WT.FirstSpike_std  = cellfun(@nanstd,analysis.WT.FirstSpike,'UniformOutput',false);
analysis.WT.FirstSpike_std  = cell2mat(cellfun(@transpose,analysis.WT.FirstSpike_std,'UniformOutput',false));

analysis.WT.NoSpikes_mean   = cellfun(@nanmean,analysis.WT.NoSpikes,'UniformOutput',false);
analysis.WT.NoSpikes_mean   = cell2mat(cellfun(@transpose,analysis.WT.NoSpikes_mean,'UniformOutput',false));

analysis.WT.NoSpikes_sem   = cellfun(@nansem,analysis.WT.NoSpikes,'UniformOutput',false);
analysis.WT.NoSpikes_sem   = cell2mat(cellfun(@transpose,analysis.WT.NoSpikes_sem,'UniformOutput',false));


analysis.WT.ISIhist_mean   = cellfun(@nanmean,analysis.WT.ISIhist,'UniformOutput',false);
analysis.WT.ISIhist_sem    = nansem(cell2mat(cellfun(@transpose,analysis.WT.ISIhist_mean(1:end),'UniformOutput',false)),2);
analysis.WT.ISIhist_mean   = nanmean(cell2mat(cellfun(@transpose,analysis.WT.ISIhist_mean(1:end),'UniformOutput',false)),2);

analysis.KO.FiringSpan = [];
analysis.KO.FirstSpike = [];
analysis.KO.NoSpikes   = [];
analysis.KO.ISIhist    = [];
for iModel = 1:noModels
    for iTrial = 1:noTrials
        temp=resmatrix.KO(iModel).firings{iFreq,iTrial};
        for iNeuron= 1:Ne
            if length(temp(temp(:,2)==iNeuron,1))>1 %lots of spikes
                analysis.KO.FiringSpan{iModel,iTrial}(iNeuron,:) = [min(temp(temp(:,2)==iNeuron,1)), max(temp(temp(:,2)==iNeuron,1))];
                analysis.KO.FirstSpike{iModel}(iTrial,iNeuron)   = min(temp(temp(:,2)==iNeuron,1));
                analysis.KO.NoSpikes{iModel}(iTrial,iNeuron)     = length(temp(temp(:,2)==iNeuron,1));
                analysis.KO.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = histc(diff(temp(temp(:,2)==iNeuron,1)),ISIbins);
                
            elseif length(temp(temp(:,2)==iNeuron,1))==1 % 1 spike
                analysis.KO.FiringSpan{iModel,iTrial}(iNeuron,:) = [min(temp(temp(:,2)==iNeuron,1)), NaN];
                analysis.KO.FirstSpike{iModel}(iTrial,iNeuron)   = min(temp(temp(:,2)==iNeuron,1));
                analysis.KO.NoSpikes{iModel}(iTrial,iNeuron)     = 1;
                analysis.KO.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = nan(size(ISIbins));
                
                
            else % no spikes
                analysis.KO.FiringSpan{iModel,iTrial}(iNeuron,:) = [NaN,NaN];
                analysis.KO.FirstSpike{iModel}(iTrial,iNeuron)   = NaN;
                analysis.KO.NoSpikes{iModel}(iTrial,iNeuron)     = 0;
                analysis.KO.ISIhist{iModel,iTrial}(iNeuron,1:length(ISIbins))    = nan(size(ISIbins));
            end
        end
    end
end
analysis.KO.FirstSpike_mean = cellfun(@nanmean,analysis.KO.FirstSpike,'UniformOutput',false);
analysis.KO.FirstSpike_mean = cell2mat(cellfun(@transpose,analysis.KO.FirstSpike_mean,'UniformOutput',false));

analysis.KO.FirstSpike_std  = cellfun(@nanstd,analysis.KO.FirstSpike,'UniformOutput',false);
analysis.KO.FirstSpike_std  = cell2mat(cellfun(@transpose,analysis.KO.FirstSpike_std,'UniformOutput',false));

analysis.KO.NoSpikes_mean   = cellfun(@nanmean,analysis.KO.NoSpikes,'UniformOutput',false);
analysis.KO.NoSpikes_mean   = cell2mat(cellfun(@transpose,analysis.KO.NoSpikes_mean,'UniformOutput',false));

analysis.KO.NoSpikes_sem   = cellfun(@nansem,analysis.KO.NoSpikes,'UniformOutput',false);
analysis.KO.NoSpikes_sem   = cell2mat(cellfun(@transpose,analysis.KO.NoSpikes_sem,'UniformOutput',false));

analysis.KO.ISIhist_mean   = cellfun(@nanmean,analysis.KO.ISIhist,'UniformOutput',false);
analysis.KO.ISIhist_sem    = nansem(cell2mat(cellfun(@transpose,analysis.KO.ISIhist_mean(1:end),'UniformOutput',false)),2);
analysis.KO.ISIhist_mean   = nanmean(cell2mat(cellfun(@transpose,analysis.KO.ISIhist_mean(1:end),'UniformOutput',false)),2);
%% spike stats plots
figure; hold on
tempWT = nanmean(histc(analysis.WT.NoSpikes_mean,0:10),2); tempWT =tempWT ./sum(tempWT);
tempKO = nanmean(histc(analysis.KO.NoSpikes_mean,0:10),2); tempKO =tempKO ./sum(tempKO);
area((0:10),tempKO,'FaceColor',[0.9 0.4 0.4],'EdgeColor',[0.1 0 0],'EdgeAlpha',0.5,'FaceAlpha',0.5,'LineWidth',1.5)
area((0:10),tempWT,'FaceColor',[0.4 0.4 0.9],'EdgeColor',[0 0 0.1],'EdgeAlpha',0.5,'FaceAlpha',0.5,'LineWidth',1.5)
xlabel('Mean no. spikes fired')
ylabel('Fraction of Ex. neurons')
%% inter-spike interval histograms
figure; hold on
plot(0:200,(histc(analysis.WT.FirstSpike_mean,0:200)),'b')
plot(0:200,(histc(analysis.KO.FirstSpike_mean,0:200)),'r')

figure; hold on
ciplot(analysis.WT.ISIhist_mean+analysis.WT.ISIhist_sem,...
    analysis.WT.ISIhist_mean-analysis.WT.ISIhist_sem,...
    ISIbins,'b')

ciplot(analysis.KO.ISIhist_mean+analysis.KO.ISIhist_sem,...
    analysis.KO.ISIhist_mean-analysis.KO.ISIhist_sem,...
    ISIbins,'r')
plot(ISIbins,analysis.WT.ISIhist_mean,'b')
plot(ISIbins,analysis.KO.ISIhist_mean,'r')
%% first spike latency vs. spike jitter
figure; hold on
iModel_WT = 3;
scatter(analysis.WT.FirstSpike_mean(1:end,iModel_WT)-10,...
    analysis.WT.NoSpikes_mean(1:end,iModel_WT),...
    analysis.WT.FirstSpike_std(1:end,iModel_WT),'b')

iModel_KO = 4;

scatter(analysis.KO.FirstSpike_mean(1:end,iModel_KO)-10,...
    analysis.KO.NoSpikes_mean(1:end,iModel_KO),...
    analysis.KO.FirstSpike_std(1:end,iModel_KO),'r')

tempWT= 0.005*(histc(analysis.WT.FirstSpike_mean(:,iModel_WT),0:1000));
tempKO= 0.005*(histc(analysis.KO.FirstSpike_mean(:,iModel_KO),0:1000));
area((0:1000)-10,tempKO,'FaceColor',[0.9 0.4 0.4],'EdgeColor',[0.1 0 0])
area((0:1000)-10,tempWT,'FaceColor',[0.4 0.4 0.9],'EdgeColor',[0 0 0.1])


fade_=[0.5 0.5 0.5];
ISI=1000*(stimFreqs(iFreq)^-1);
scatter([0, ISI, 2*ISI, 3*ISI, 4*ISI], 0.8*[1 1 1 1 1], 30*[1 1 1 1 1], repmat(fade_,5,1),'^','LineWidth',1.5)
text(2*ISI,0.5,'3rd','HorizontalAlignment','center','Color',fade_)
text(3*ISI,0.5,'4th','HorizontalAlignment','center','Color',fade_)
text(4*ISI,0.5,'5th','HorizontalAlignment','center','Color',fade_)
text(4.1*ISI,0.5,'stimulus','HorizontalAlignment','left','Color',fade_)

text(87,5.6,'First Spike Jitter','HorizontalAlignment','center','Color',fade_)
scatter(82,4.5,5,[0.6 0.6 0.6]);    text(92,4.52,'5ms','Color',fade_,'HorizontalAlignment','right')
scatter(82,4.73,10,[0.6 0.6 0.6]);  text(92,4.77,'10ms','Color',fade_,'HorizontalAlignment','right')
scatter(82,4.95,50,[0.6 0.6 0.6]);     text(92,5.01,'50ms','Color',fade_,'HorizontalAlignment','right')
scatter(82,5.25,100,[0.6 0.6 0.6]); text(92,5.27,'100ms','Color',fade_,'HorizontalAlignment','right')
% rectangle('Position',[78,4.3,18,1.5],'EdgeColor',[0.6 0.6 0.6])
ylabel('Av. no. spikes fired')
xlabel('Mean first spike latency (ms)')
axis([25 100 0 6 ])
