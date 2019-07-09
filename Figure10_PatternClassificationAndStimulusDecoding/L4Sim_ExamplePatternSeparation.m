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
% This script recapitulates the simulation analysis results for Fig 10A-C in the above paper.
% This code plots example spike trains from single units in response to inputs with/without added distractors.
% Quantification of dpike train modulation (change in rate, first spike
% time) is provided, as well as through a spike train metric 
% (van Rossum MC, Neural computation 2001 Apr;13(4):751-63)
% 
% This code uses an optimised implementation of the van Rossum distance
% (Haughton C, Kreuz T, Network. 2012;23(1-2):48-58.). This dependency
% function is available online:
% provided (MATLAB) by Thomas Kreuz at: http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/VanRossum.html
% and (C++) by Conor Haughton at: https://sourceforge.net/p/spikemetrics/code/ci/default/tree/
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

%% Setup
clear; close all

if ispc
    pat= 'E:\ModellingResults\PatternSeparation\';
else ismac
    pat = '/Volumes/5TB external/ModellingResults/PatternSeparation/';
end

stimFreqs=[5 10 20 50];
stimPosition = [0 1 2 3 4];
noModels=5;
noStimConds=length(stimPosition);
noTrials=10;
Ne=800;
%% Import WT
cd(pat)
A  = dir([pat 'WT_20Hz*']);
cd(A(1).name)

for iTrial =1:noTrials
    iTrial
    B(1) = load(sprintf('WT_addedDistractor_5_0_%d.mat',iTrial)); % model/stim/trial
    B(2) = load(sprintf('WT_addedDistractor_5_2_%d.mat',iTrial)); % model/stim/trial
    WTno_spikes{1}(iTrial,:)= sum(B(1).results.spike_raster(1:Ne,:),2);
    WTno_spikes{2}(iTrial,:)= sum(B(2).results.spike_raster(1:Ne,:),2);
end
%% Import KO
cd(pat)
A  = dir([pat 'KO_20Hz*']);
cd(A(1).name)

for iTrial =1:noTrials
    iTrial
    C(1) = load(sprintf('KO_addedDistractor_2_0_%d.mat',iTrial)); % model/stim/trial
    C(2) = load(sprintf('KO_addedDistractor_2_2_%d.mat',iTrial)); % model/stim/trial
    KOno_spikes{1}(iTrial,:)= sum(C(1).results.spike_raster(1:Ne,:),2);
    KOno_spikes{2}(iTrial,:)= sum(C(2).results.spike_raster(1:Ne,:),2);
end
%% Mean change in spike no.
figure; hold on
XY_max = 7;
subplot(2,1,1); hold on
plot([0 XY_max],[0 XY_max],':k')
axis([0 XY_max 0 XY_max ])
axis square
scatter(mean(WTno_spikes{1}),mean(WTno_spikes{2}),10,'b','filled')
subplot(2,1,2); hold on
scatter(mean(KOno_spikes{1}),mean(KOno_spikes{2}),10,'r','filled')
plot([0 XY_max],[0 XY_max],':k')
axis([0 XY_max 0 XY_max ])
axis square
% xlabel('Pattern 1 response (no. spikes)')
% ylabel('Pattern 2 response (no. spikes)')
%% Mean change in spike no.
figure; hold on
XY_max = 7;
subplot(2,1,1); hold on
plot([0 XY_max],[0 XY_max],':k')
axis([0 XY_max 0 XY_max ])
axis square
scatter(WTno_spikes{1}(:),WTno_spikes{2}(:),20,'b','filled')
subplot(2,1,2); hold on
scatter(KOno_spikes{1}(:),KOno_spikes{2}(:)+0.1,20,'r','filled')
plot([0 XY_max],[0 XY_max],':k')
axis([0 XY_max 0 XY_max ])
axis square
% xlabel('Pattern 1 response (no. spikes)')
% ylabel('Pattern 2 response (no. spikes)')
%% 2D histogram of change in spike no.
counts = 0:10;
for i=0:10
    for j=0:10
        WTno_spikesCollapsed(i+1,j+1) = sum(WTno_spikes{1}(:)==i & WTno_spikes{2}(:) ==j);
        KOno_spikesCollapsed(i+1,j+1) = sum(KOno_spikes{1}(:)==i & KOno_spikes{2}(:) ==j);
        xvals(i+1,j+1) = i+1;
        yvals(i+1,j+1) = j+1;
    end
end
WTno_spikesCollapsed=WTno_spikesCollapsed./Ne;
KOno_spikesCollapsed=KOno_spikesCollapsed./Ne;
WTno_spikesCollapsed(WTno_spikesCollapsed==0)=NaN;
KOno_spikesCollapsed(KOno_spikesCollapsed==0)=NaN;
figure;
subplot(1,2,1); hold on
imagesc(counts,counts,WTno_spikesCollapsed'); set(gca,'Ydir','normal')
plot([0 10],[0 10],':k')
axis([0 10 0 10])
subplot(1,2,2); hold on
imagesc(counts,counts,KOno_spikesCollapsed'); set(gca,'Ydir','normal')
plot([0 10],[0 10],':k')
colormap([1 1 1;jet])
axis([0 10 0 10])

figure;
subplot(1,2,1); hold on
scatter(xvals(:),yvals(:),WTno_spikesCollapsed(:)*100,'b','filled','MarkerEdgeColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
plot([0 10],[0 10],':k')
axis([0 10 0 10])
axis square
subplot(1,2,2); hold on
scatter(xvals(:),yvals(:),KOno_spikesCollapsed(:)*100,'r','filled','MarkerEdgeColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
plot([0 10],[0 10],':k')
colormap([1 1 1;jet])
axis([0 10 0 10])
axis square
%% Get pairwise spiketrain similarity using van Rossum distance
Ne=800;
for i =1:Ne
    
    X = B(1).results.spike_raster(i,:);
    Y = B(2).results.spike_raster(i,:);
    Disc_WT(i)= abs(vanRossum(find(X)/2,find(Y)/2,50));
    
    
    X = C(1).results.spike_raster(i,:);
    Y = C(2).results.spike_raster(i,:);
    Disc_KO(i)= abs(vanRossum(find(X)/2,find(Y)/2,50));
end
bins = 0:0.01:1;
figure; hold on
y=smooth_hist(histc(sort(Disc_WT),bins));
stairs(bins,y./sum(y)*100,'b','LineWidth',1.5)
y=smooth_hist(histc(sort(Disc_KO),bins));
stairs(bins,y./sum(y)*100,'r','LineWidth',1.5)
%% shift in no spikes/first spike time vs pairwise spike train similarity

figure; hold on
X = B(1).results.time_stats.spike_times(1:Ne);
idx = ~cellfun('isempty',X);
firstspike_norm = nan(size(X));
firstspike_norm(idx) = cellfun(@(x)x(1),X(idx));
no_spikes_norm = cellfun(@numel,X);
X = B(2).results.time_stats.spike_times(1:Ne);
idx = ~cellfun('isempty',X);
firstspike_dist = nan(size(X));
firstspike_dist(idx) = cellfun(@(x)x(1),X(idx));
no_spikes_dist = cellfun(@numel,X);
delta_firstspike = firstspike_dist-firstspike_norm;
delta_noSpikes = no_spikes_dist-no_spikes_norm;

scatter(delta_firstspike,delta_noSpikes,1+(20*Disc_WT).^2,'b','filled','MarkerEdgeColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)

X = C(1).results.time_stats.spike_times(1:Ne);
idx = ~cellfun('isempty',X);
firstspike_norm = nan(size(X));
firstspike_norm(idx) = cellfun(@(x)x(1),X(idx));
no_spikes_norm = cellfun(@numel,X);
X = C(2).results.time_stats.spike_times(1:Ne);
idx = ~cellfun('isempty',X);
firstspike_dist = nan(size(X));
firstspike_dist(idx) = cellfun(@(x)x(1),X(idx));
no_spikes_dist = cellfun(@numel,X);
delta_firstspike = firstspike_dist-firstspike_norm;
delta_noSpikes = no_spikes_dist-no_spikes_norm;

scatter(delta_firstspike,delta_noSpikes,1+(20*Disc_KO).^2,'r','filled','MarkerEdgeColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
ylims = get(gca,'YLim');xlims = get(gca,'XLim');
plot([0 0],ylims,':k')
plot(xlims,[0 0],':k')
colormap hot
%% Example WT membrane potential
tb = B(1).M.p.tb;
Ne = B(1).M.p.Ne;
figure
subplot(1,2,1); hold on

% lowestD = randsample(Ne,1);
% lowestD = find(Disc_WT==min(Disc_WT));
% lowestD = FindClosestIndex(Disc_WT,0.57);
lowestD = find(sum(B(1).results.spike_raster((1:Ne),:),2)==1 & ...
    sum(B(2).results.spike_raster((1:Ne),:),2)==2);
lowestD = lowestD(2);
plot(tb,B(2).results.V_out(lowestD,:),'color',[0.9 0.6 0.3],'LineWidth',1.5)
plot(tb,B(1).results.V_out(lowestD,:),'color','k','LineWidth',1.5)

stim_times = B(2).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color',[0.9 0.6 0.3],'LineWidth',1.5)
stim_times = B(1).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color','k','LineWidth',1.5)
axis([-10 500 -90 25])

subplot(1,2,2); hold on
% lowestD = randsample(Ne,1);
% lowestD = find(Disc_WT==min(Disc_WT));
% lowestD = FindClosestIndex(Disc_WT,0.57);
lowestD = find(sum(B(1).results.spike_raster((1:Ne),:),2)==2 & ...
    sum(B(2).results.spike_raster((1:Ne),:),2)==3);
lowestD = lowestD(1);
plot(tb,B(2).results.V_out(lowestD,:),'color',[0.9 0.6 0.3],'LineWidth',1.5)
plot(tb,B(1).results.V_out(lowestD,:),'color','k','LineWidth',1.5)

stim_times = B(2).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color',[0.9 0.6 0.3],'LineWidth',1.5)
stim_times = B(1).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color','k','LineWidth',1.5)
axis([-10 500 -90 25])
%% Example KO membrane potential
tb = B(1).M.p.tb;
Ne = B(1).M.p.Ne;
figure
subplot(1,2,1); hold on

% lowestD = randsample(Ne,1);
% lowestD = find(Disc_WT==min(Disc_WT));
% lowestD = FindClosestIndex(Disc_WT,0.57);
lowestD = find(sum(C(1).results.spike_raster((1:Ne),:),2)==3 & ...
    sum(C(2).results.spike_raster((1:Ne),:),2)==3);
lowestD = lowestD(1);
plot(tb,C(2).results.V_out(lowestD,:),'color',[0.9 0.6 0.3],'LineWidth',1.5)
plot(tb,C(1).results.V_out(lowestD,:),'color','k','LineWidth',1.5)

stim_times = C(2).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color',[0.9 0.6 0.3],'LineWidth',1.5)
stim_times = C(1).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color','k','LineWidth',1.5)
axis([-10 500 -90 25])

subplot(1,2,2); hold on
% lowestD = randsample(Ne,1);
% lowestD = find(Disc_WT==min(Disc_WT));
% lowestD = FindClosestIndex(Disc_WT,0.57);
lowestD = find(sum(C(1).results.spike_raster((1:Ne),:),2)==6 & ...
    sum(C(2).results.spike_raster((1:Ne),:),2)==6);
lowestD = lowestD(1);
plot(tb,C(2).results.V_out(lowestD,:),'color',[0.9 0.6 0.3],'LineWidth',1.5)
plot(tb,C(1).results.V_out(lowestD,:),'color','k','LineWidth',1.5)

stim_times = C(2).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color',[0.9 0.6 0.3],'LineWidth',1.5)
stim_times = C(1).M.input.ZAP.stim_times;
plot([stim_times;stim_times],[-75*ones(1,length(stim_times));-70*ones(1,length(stim_times))],'color','k','LineWidth',1.5)
axis([-10 500 -90 25])
%% Load all trials (WT)

cd(pat)
A  = dir([pat 'WT_20Hz*']);
cd(A(1).name)

for iTrial =1:noTrials
    iTrial
    B(1) = load(sprintf('WT_addedDistractor_5_0_%d.mat',iTrial)); % model/stim/trial
    B(2) = load(sprintf('WT_addedDistractor_5_2_%d.mat',iTrial)); % model/stim/trial
    
    X = B(1).results.time_stats.spike_times(1:Ne);
    idx = ~cellfun('isempty',X);
    firstspike_norm = nan(size(X));
    firstspike_norm(idx) = cellfun(@(x)x(1),X(idx));
    no_spikes_norm = cellfun(@numel,X);
    X = B(2).results.time_stats.spike_times(1:Ne);
    idx = ~cellfun('isempty',X);
    firstspike_dist = nan(size(X));
    firstspike_dist(idx) = cellfun(@(x)x(1),X(idx));
    no_spikes_dist = cellfun(@numel,X);
    delta_firstspikeWT(:,iTrial) = firstspike_dist-firstspike_norm;
    delta_noSpikesWT(:,iTrial) = no_spikes_dist-no_spikes_norm;
    
    
    for i =1:Ne
        
        X = B(1).results.spike_raster(i,:);
        Y = B(2).results.spike_raster(i,:);
        Disc_WT(i,iTrial)= abs(vanRossum(find(X)/2,find(Y)/2,50));
        
    end
end
%% Load all trials (KO)


cd(pat)
A  = dir([pat 'KO_20Hz*']);
cd(A(1).name)

for iTrial =1:noTrials
    iTrial
    C(1) = load(sprintf('KO_addedDistractor_2_0_%d.mat',iTrial)); % model/stim/trial
    C(2) = load(sprintf('KO_addedDistractor_2_2_%d.mat',iTrial)); % model/stim/trial
    
    X = C(1).results.time_stats.spike_times(1:Ne);
    idx = ~cellfun('isempty',X);
    firstspike_norm = nan(size(X));
    firstspike_norm(idx) = cellfun(@(x)x(1),X(idx));
    no_spikes_norm = cellfun(@numel,X);
    X = C(2).results.time_stats.spike_times(1:Ne);
    idx = ~cellfun('isempty',X);
    firstspike_dist = nan(size(X));
    firstspike_dist(idx) = cellfun(@(x)x(1),X(idx));
    no_spikes_dist = cellfun(@numel,X);
    delta_firstspikeKO(:,iTrial) = firstspike_dist-firstspike_norm;
    delta_noSpikesKO(:,iTrial) = no_spikes_dist-no_spikes_norm;
    
    
    for i =1:Ne
        
        X = C(1).results.spike_raster(i,:);
        Y = C(2).results.spike_raster(i,:);
        Disc_KO(i,iTrial)= abs(vanRossum(find(X)/2,find(Y)/2,50));
        
    end
end
%% Shift in no spikes/first spike time vs pairwise spike train similarity - all trials

figure; hold on

scatter(nanmean(delta_firstspikeWT,2),nanmean(delta_noSpikesWT,2),1+100*(20*nanmean(Disc_WT,2)).^2,'b','filled','MarkerEdgeColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
scatter(nanmean(delta_firstspikeKO,2),nanmean(delta_noSpikesKO,2),1+(20*nanmean(Disc_KO,2)).^2,'r','filled','MarkerEdgeColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
% scatter(delta_firstspikeWT(:,10),delta_noSpikesWT(:,10),1+(20*Disc_WT(:,10)).^2,'b','filled','MarkerEdgeColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
% scatter(nanmean(delta_firstspikeKO,2),nanmean(delta_noSpikesKO,2),1+(20*nanmean(Disc_KO,2)).^2,'r','filled','MarkerEdgeColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.8)
ylims = get(gca,'YLim');xlims = get(gca,'XLim');
plot([0 0],ylims,':k')
plot(xlims,[0 0],':k')
colormap hot
%% Histgrams of distractor discriminability by genotype

bins = 0:0.01:1;
figure; hold on
y=smooth_hist(histc(mean(Disc_WT,2),bins));
stairs(bins,y./sum(y)*100,'b','LineWidth',1.5)
y=smooth_hist(histc(mean(Disc_KO,2),bins));
stairs(bins,y./sum(y)*100,'r','LineWidth',1.5)
xlabel({'Spike train Dissimilarity','(van Rossum distance)'})
ylabel('% of neurons')