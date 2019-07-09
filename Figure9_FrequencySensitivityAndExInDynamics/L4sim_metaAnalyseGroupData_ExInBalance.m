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
% Fig 9E-F in the above paper and plots example firing rates of Ex and In 
% pools of neurons and relative timing of group firing.
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
%% set parameters
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

%% Import group data
cd(recFolderWT)
for iModel=1:noModels%noModels;  
    for iFreq = 4%1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('WT model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            meta.dynamics.SpikeDensity.mean_PDF_ex{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_ex;
            meta.dynamics.SpikeDensity.mean_PDF_in{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_in;
            meta.time_stats.hist_bins{iModel}                           = results.time_stats.hist_bins;
            meta.time_stats.spike_rate_mean_Ex{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_Ex;
            meta.time_stats.spike_rate_mean_In{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_In;
            temp=results.firings;
            meta.firings_E{iModel}{trial_id}=results.firings_E;
            meta.firings_I{iModel}{trial_id}=results.firings_I;
            temp2=results.time_stats.spike_rate;
            for n_id=1:M.p.Ne+M.p.Ni
                spike_times{n_id,trial_id}=temp(temp(:,2)==n_id,1);
                spike_rate(n_id,trial_id)              = temp2(n_id);
            end
            M.p.t_max=1000;
            % Spike density functions
            SDF_temp = SDF(results.spike_raster,0.005,[1 M.p.t_max],10);
            mean_SDF_temp(:,trial_id)    = nanmean(SDF_temp.PDF_trimmed(1:M.p.Ne,:),1);
        end
    end
    
    meta.SDF.mean_SDF_Ex{iModel}            = mean_SDF_temp;
    meta.SDF.grand_mean_SDF_Ex(iModel,:)    = nanmean(mean_SDF_temp,2);
    meta.SDF.grand_sem_SDF_Ex(iModel,:)     = nansem(mean_SDF_temp,2);
    meta.time_stats.spike_rate{iModel}      = spike_rate';
meta.time_stats.spike_rate_mean(iModel,:) = nanmean(meta.time_stats.spike_rate{iModel},1); 
% meta.time_stats.spike_rate_mean(isnan(meta.time_stats.spike_rate_mean))=0.1;
meta.connectivity{iModel}               = results.connectivity;
meta.spike_times{iModel}                = spike_times;


    for iFreq = 4%1:noStimConds
        for n_id=1:M.p.Ne
            meta.times.count.mean(iModel,n_id)  = nanmean(cellfun(@numel,meta.spike_times{iModel}( n_id,:)));
            meta.times.count.CoV(iModel,n_id)   = nanstd(cellfun(@numel,meta.spike_times{iModel}(n_id,:)))./nanmean(cellfun(@numel,meta.spike_times{iModel}(n_id,:)));
            for trial_id=1:noTrials
                temp=meta.spike_times{iModel}{n_id,trial_id};
                try
                    meta.times.first{iModel}(n_id,trial_id)=temp(1);
                catch
                    meta.times.first{iModel}(n_id,trial_id)=NaN;
                end
            end
        end
    end
    meta.times.first_mean{iModel} = nanmean(meta.times.first{iModel},2);
    meta.times.first_std{iModel}  = nanstd(meta.times.first{iModel},2);
    meta.times.first_CoV{iModel}  = meta.times.first_std{iModel}./meta.times.first_mean{iModel};
end
WTmeta=meta;
cd(recFolderKO)
for iModel=1:noModels%noModels;  
    for iFreq = 4%1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('KO model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            meta.dynamics.SpikeDensity.mean_PDF_ex{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_ex;
            meta.dynamics.SpikeDensity.mean_PDF_in{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_in;
            meta.time_stats.hist_bins{iModel}                           = results.time_stats.hist_bins;
            meta.time_stats.spike_rate_mean_Ex{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_Ex;
            meta.time_stats.spike_rate_mean_In{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_In;
            temp=results.firings;
            meta.firings_E{iModel}{trial_id}=results.firings_E;
            meta.firings_I{iModel}{trial_id}=results.firings_I;
            temp2=results.time_stats.spike_rate;
            for n_id=1:M.p.Ne+M.p.Ni
                spike_times{n_id,trial_id}=temp(temp(:,2)==n_id,1);
                spike_rate(n_id,trial_id)              = temp2(n_id);
            end
            M.p.t_max=1000;
            % Spike density functions
            SDF_temp = SDF(results.spike_raster,0.005,[1 M.p.t_max],10);
            mean_SDF_temp(:,trial_id)    = nanmean(SDF_temp.PDF_trimmed(1:M.p.Ne,:),1);
        end
    end
    
    meta.SDF.mean_SDF_Ex{iModel}            = mean_SDF_temp;
    meta.SDF.grand_mean_SDF_Ex(iModel,:)    = nanmean(mean_SDF_temp,2);
    meta.SDF.grand_sem_SDF_Ex(iModel,:)     = nansem(mean_SDF_temp,2);
    meta.time_stats.spike_rate{iModel}      = spike_rate';
meta.time_stats.spike_rate_mean(iModel,:) = nanmean(meta.time_stats.spike_rate{iModel},1); 
meta.connectivity{iModel}               = results.connectivity;
meta.spike_times{iModel}                = spike_times;


    for iFreq = 4%1:noStimConds
        for n_id=1:M.p.Ne
            meta.times.count.mean(iModel,n_id)  = nanmean(cellfun(@numel,meta.spike_times{iModel}( n_id,:)));
            meta.times.count.CoV(iModel,n_id)   = nanstd(cellfun(@numel,meta.spike_times{iModel}(n_id,:)))./nanmean(cellfun(@numel,meta.spike_times{iModel}(n_id,:)));
            for trial_id=1:noTrials
                temp=meta.spike_times{iModel}{n_id,trial_id};
                try
                    meta.times.first{iModel}(n_id,trial_id)=temp(1);
                catch
                    meta.times.first{iModel}(n_id,trial_id)=NaN;

cd(recFolderKO)
for iModel=1:noModels%noModels;  
    for iFreq = 4%1:noStimConds
        for iTrial = 1:noTrials
            disp(sprintf('KO model no. %d, stim frequency: %dHz, trial no. %d',iModel,stimFreqs(iFreq),iTrial))
            A = dir(['*' num2str(iModel) '_' num2str(stimFreqs(iFreq)) '_' num2str(iTrial) '.mat']);
            load(A(1).name);
            meta.EIcorr.lags                                            = results.EIcorr.lags;
            meta.EIcorr.r{iModel}(:,trial_id)                           = results.EIcorr.r;
            meta.spectro.FFT.freq{iModel}                               = results.spectro.FFT.freq;
            meta.spectro.FFT.power{iModel}(:,trial_id)                  = results.spectro.FFT.power;
            meta.dynamics.SpikeDensity.mean_PDF_ex{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_ex;
            meta.dynamics.SpikeDensity.mean_PDF_in{iModel}(:,trial_id)  = results.dynamics.SpikeDensity.mean_PDF_in;
            meta.time_stats.hist_bins{iModel}                           = results.time_stats.hist_bins;
            meta.time_stats.spike_rate_mean_Ex{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_Ex;
            meta.time_stats.spike_rate_mean_In{iModel}(:,trial_id)      = results.time_stats.spike_rate_mean_In;
            meta.mean_Ex_to_Ex{iModel}(:,trial_id)                      = results.mean_Ex_to_Ex;
            meta.mean_In_to_Ex{iModel}(:,trial_id)                      = results.mean_In_to_Ex;
            temp=results.firings;
            meta.firings_E{iModel}{trial_id}=results.firings_E;
            meta.firings_I{iModel}{trial_id}=results.firings_I;
            temp2=results.time_stats.spike_rate;
            for n_id=1:M.p.Ne+M.p.Ni
                spike_times{n_id,trial_id}=temp(temp(:,2)==n_id,1);
                spike_rate(n_id,trial_id)              = temp2(n_id);
            end
            M.p.t_max=1000;
            % Spike density functions
            SDF_temp = SDF(results.spike_raster,0.005,[1 M.p.t_max],10);
            mean_SDF_temp(:,trial_id)    = nanmean(SDF_temp.PDF_trimmed(1:M.p.Ne,:),1);
        end
    end
    
    meta.SDF.mean_SDF_Ex{iModel}            = mean_SDF_temp;
    meta.SDF.grand_mean_SDF_Ex(iModel,:)    = nanmean(mean_SDF_temp,2);
    meta.SDF.grand_sem_SDF_Ex(iModel,:)     = nansem(mean_SDF_temp,2);
    meta.time_stats.spike_rate{iModel}      = spike_rate';
meta.time_stats.spike_rate_mean(iModel,:) = nanmean(meta.time_stats.spike_rate{iModel},1); 
meta.connectivity{iModel}               = results.connectivity;
meta.spike_times{iModel}                = spike_times;


    for iFreq = 4%1:noStimConds
        for n_id=1:M.p.Ne
            meta.times.count.mean(iModel,n_id)  = nanmean(cellfun(@numel,meta.spike_times{iModel}( n_id,:)));
            meta.times.count.CoV(iModel,n_id)   = nanstd(cellfun(@numel,meta.spike_times{iModel}(n_id,:)))./nanmean(cellfun(@numel,meta.spike_times{iModel}(n_id,:)));
            for trial_id=1:noTrials
                temp=meta.spike_times{iModel}{n_id,trial_id};
                try
                    meta.times.first{iModel}(n_id,trial_id)=temp(1);
                catch
                    meta.times.first{iModel}(n_id,trial_id)=NaN;
                end
            end
        end
    end
    meta.times.first_mean{iModel} = nanmean(meta.times.first{iModel},2);
    meta.times.first_std{iModel}  = nanstd(meta.times.first{iModel},2);
    meta.times.first_CoV{iModel}  = meta.times.first_std{iModel}./meta.times.first_mean{iModel};
    meta.spectro.FFT.power_mean=nanmean(meta.spectro.FFT.power{iModel},2);
    meta.spectro.FFT.power_SEM=nansem(meta.spectro.FFT.power{iModel},2);
end
                end
            end
        end
    end
    meta.times.first_mean{iModel} = nanmean(meta.times.first{iModel},2);
    meta.times.first_std{iModel}  = nanstd(meta.times.first{iModel},2);
    meta.times.first_CoV{iModel}  = meta.times.first_std{iModel}./meta.times.first_mean{iModel};

end
KOmeta=meta;
clear trial_id results Fname temp temp2 n_id spike_times spike_rate meta
%% Collapse across models

WTmeta.dynamics.SpikeDensity.mean_PDF_ex_=cell2mat(WTmeta.dynamics.SpikeDensity.mean_PDF_ex);
WTmeta.dynamics.SpikeDensity.mean_PDF_ex_mean = nanmean(WTmeta.dynamics.SpikeDensity.mean_PDF_ex_,2);
WTmeta.dynamics.SpikeDensity.mean_PDF_ex_sem = nansem(WTmeta.dynamics.SpikeDensity.mean_PDF_ex_,2);

KOmeta.dynamics.SpikeDensity.mean_PDF_ex_=cell2mat(KOmeta.dynamics.SpikeDensity.mean_PDF_ex);
KOmeta.dynamics.SpikeDensity.mean_PDF_ex_mean = nanmean(KOmeta.dynamics.SpikeDensity.mean_PDF_ex_,2);
KOmeta.dynamics.SpikeDensity.mean_PDF_ex_sem = nansem(KOmeta.dynamics.SpikeDensity.mean_PDF_ex_,2);

WTmeta.dynamics.SpikeDensity.mean_PDF_in_=cell2mat(WTmeta.dynamics.SpikeDensity.mean_PDF_in);
WTmeta.dynamics.SpikeDensity.mean_PDF_in_mean = nanmean(WTmeta.dynamics.SpikeDensity.mean_PDF_in_,2);
WTmeta.dynamics.SpikeDensity.mean_PDF_in_sem = nansem(WTmeta.dynamics.SpikeDensity.mean_PDF_in_,2);

KOmeta.dynamics.SpikeDensity.mean_PDF_in_=cell2mat(KOmeta.dynamics.SpikeDensity.mean_PDF_in);
KOmeta.dynamics.SpikeDensity.mean_PDF_in_mean = nanmean(KOmeta.dynamics.SpikeDensity.mean_PDF_in_,2);
KOmeta.dynamics.SpikeDensity.mean_PDF_in_sem = nansem(KOmeta.dynamics.SpikeDensity.mean_PDF_in_,2);

WTmeta.mean_Ex_to_Ex_=cell2mat(WTmeta.mean_Ex_to_Ex);
WTmeta.mean_Ex_to_Ex_mean=nanmean(WTmeta.mean_Ex_to_Ex_,2);
WTmeta.mean_Ex_to_Ex_sem=nansem(WTmeta.mean_Ex_to_Ex_,2);
KOmeta.mean_Ex_to_Ex_=cell2mat(KOmeta.mean_Ex_to_Ex);
KOmeta.mean_Ex_to_Ex_mean=nanmean(KOmeta.mean_Ex_to_Ex_,2);
KOmeta.mean_Ex_to_Ex_sem=nansem(KOmeta.mean_Ex_to_Ex_,2);

WTmeta.mean_In_to_Ex_=cell2mat(WTmeta.mean_In_to_Ex);
WTmeta.mean_In_to_Ex_mean=nanmean(WTmeta.mean_In_to_Ex_,2);
WTmeta.mean_In_to_Ex_sem=nansem(WTmeta.mean_In_to_Ex_,2);
KOmeta.mean_In_to_Ex_=cell2mat(KOmeta.mean_In_to_Ex);
KOmeta.mean_In_to_Ex_mean=nanmean(KOmeta.mean_In_to_Ex_,2);
KOmeta.mean_In_to_Ex_sem=nansem(KOmeta.mean_In_to_Ex_,2);

WTmeta.times.first_CoV_=cell2mat(WTmeta.times.first_CoV);
WTmeta.times.first_CoV_=WTmeta.times.first_CoV_(1:end);
WTmeta.times.count_=WTmeta.times.count.mean(1:end);
KOmeta.times.first_CoV_=cell2mat(KOmeta.times.first_CoV);
KOmeta.times.count_=KOmeta.times.count.mean(1:end);
KOmeta.times.first_CoV_=KOmeta.times.first_CoV_(1:end);
%% Plot EI SDFs
dt=1;
t_max=length(WTmeta.dynamics.SpikeDensity.mean_PDF_ex_mean);
ISI=1000*(stimFreqs(iFreq)^-1);
% figure; hold on
%     plot((1:t_max)*dt,100*WTmeta.dynamics.SpikeDensity.mean_PDF_ex_mean,'b','LineWidth',1.5);
%     plot((1:t_max)*dt,100*KOmeta.dynamics.SpikeDensity.mean_PDF_ex_mean,'r','LineWidth',1.5);
%     axis([0 500 0 Inf])
figure; 
subplot(2,1,1);hold on
    plot((1:t_max)*dt,100*WTmeta.dynamics.SpikeDensity.mean_PDF_ex_mean,'color',[0.6 0.6 0.6],'LineWidth',2);
    plot((1:t_max)*dt,100*WTmeta.dynamics.SpikeDensity.mean_PDF_in_mean,'color',[0 1 0],'LineWidth',2);
                    scatter([10,...
                      10+ISI,...
                      10+2*ISI,...
                      10+3*ISI,...
                      10+4*ISI], [0 0 0 0 0],'.k')
                  axis([0 500 0 25])
                  axis off
legend('Excitatory neurons','Inhibitory neurons'); legend('boxoff')
subplot(2,1,2);hold on
    plot((1:t_max)*dt,100*KOmeta.dynamics.SpikeDensity.mean_PDF_ex_mean,'color',[0.6 0.6 0.6],'LineWidth',2);
    plot((1:t_max)*dt,100*KOmeta.dynamics.SpikeDensity.mean_PDF_in_mean,'color',[0 1 0],'LineWidth',2);
    scatter([10,...
                      10+ISI,...
                      10+2*ISI,...
                      10+3*ISI,...
                      10+4*ISI], [0 0 0 0 0],'.k')
                  
                  
plot([250 300],[5 5],'k','LineWidth',1.5)
plot([250 250],[5 10],'k','LineWidth',1.5)
text(380, 7,{'5% synchronous';'50ms'})

    axis([0 500 0 25])
    axis off
%% Plot Ex/In synchrony
figure; hold on
twin=50:1000;
plot(100*WTmeta.dynamics.SpikeDensity.mean_PDF_ex_mean(twin),100*WTmeta.dynamics.SpikeDensity.mean_PDF_in_mean(twin),'b','LineWidth',2);
plot(100*KOmeta.dynamics.SpikeDensity.mean_PDF_ex_mean(twin),100*KOmeta.dynamics.SpikeDensity.mean_PDF_in_mean(twin),'r','LineWidth',2);
plot([0 25],[0 25],':k')
axis([0 25 0 25])
% title ('Synchronous network fraction')
xlabel('Synchronous Ex. neurons (%)')
ylabel('Synchronous In. neurons (%)')
axis square
