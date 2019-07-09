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
% - This loop simulates multiple repetitive trials of model TC stimulation for one random network seed.
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


%% Outer loop: select genotype to use model parameters for:
clear;close all      
Genotype = 'WT'; % {'WT','KO'};
%% setup network and synaptic parameters

switch Genotype
    case 'WT'
        
        
    M    = L4sim_DesignNetwork(800,150,...                 % number of Ex, In neurons
                                0.7,2000,...               % TC->In STP parameters: dep. per stim and recovery time constant
                                0.7,500,...                % TC->Ex STP parameters: dep. per stim and recovery time constant
                                0.6,500,...                % In->Ex STP parameters: dep. per stim and recovery time constant
                                0.45,500,...               % Ex->In STP parameters: dep. per stim and recovery time constant
                                0.6,500,...                % In->In STP parameters: dep. per stim and recovery time constant
                                0.6,2000,...               % Ex->Ex STP parameters: dep. per stim and recovery time constant
                                [0.18,0.0001,0.0001],...   % [p(Ex>Ex), mean g_max AMPA, mean g_max NMDA]           
                                [0.5, 0.0001],...          % [p(Ex>In), mean strength]           
                                [0.55,0.0002],...          % [p(In>Ex), mean strength]           
                                [0.5, 0.0005],...          % [p(In>In), mean strength]  
                                [0.005,0.025,0.005],...    % external gmax: AMPA[Ex, In], NMDA[Ex] (nS)
                                1,'WT');                   % plot online, genotype ('WT' or 'KO')   
                        
    case 'KO'
    M    = Iz_design_con_matrix_Fmr1_COBN(800,150,...  
                                0.2,8000,...           
                                0.5,500,...            
                                0.6,3000,...           
                                0.2,500,...            
                                0.6,3000,...           
                                0.3,5e3,...            
                                [0.18,0.0001,0.0001],...
                                [0.2,0.0001],...        
                                [0.3,0.0002],...        
                                [0.3,0.0005],...        
                                [0.005,0.025,0.005],... 
                                1,'KO');  
end
%% Design stimulation parameters
M    = L4sim_MakePulseInput(M,200,0.5,[1 1],[5 50], 1); % (M,t_max,dt,Amp,StimPattern,plot_online)
%% run simulation loop
no_trials=10;
for iTrial=1:no_trials
    
    % randomly seed Vm to introduce trial-trial variability
    M.p.vr = M.p.vr + 0.5*randn(M.p.Ne+M.p.Ni,1); 
    r=randperm(M.p.Ne)'; % Shuffle cells receiving TC stim
    
    M.net.S_ext(1:M.p.Ne)=M.net.S_ext(r);
    M.net.S_ext_NMDA(1:M.p.Ne)=M.net.S_ext_NMDA(r);
    
    results  = L4sim_RunModel(M,0,1); % (Model, online plot type, turn on shot-term plasticity Y/N)
    results  = L4sim_Analyse(M,results);
    L4sim_PlotResults(M,results);
    Fname = sprintf('%s_trial_%d.mat',Genotype,iTrial);

    save(Fname)
    clear results
    close all
end; clear iTrial
save(sprintf('%s_meta.mat',Genotype))
%% meta-analysis
% no_trials=50;
for iTrial=1:no_trials
    iTrial
    load(sprintf('%s_trial_%d.mat',Genotype,iTrial));    
    meta.EIcorr.lags                                    = results.EIcorr.lags;
    meta.EIcorr.r(:,iTrial)                           = results.EIcorr.r;
    meta.spectro.FFT.freq                               = results.spectro.FFT.freq;
    meta.spectro.FFT.power(:,iTrial)                  = results.spectro.FFT.power;
    meta.dynamics.SpikeDensity.mean_PDF_ex(:,iTrial)  = results.dynamics.SpikeDensity.mean_PDF_ex;
    meta.dynamics.SpikeDensity.mean_PDF_in(:,iTrial)  = results.dynamics.SpikeDensity.mean_PDF_in;
    meta.time_stats.hist_bins                           = results.time_stats.hist_bins;
    meta.time_stats.spike_rate_mean_Ex(:,iTrial)      = results.time_stats.spike_rate_mean_Ex;
    meta.time_stats.spike_rate_mean_In(:,iTrial)      = results.time_stats.spike_rate_mean_In;
    meta.mean_Ex_to_Ex(:,iTrial)                      = results.mean_Ex_to_Ex;
    meta.mean_In_to_Ex(:,iTrial)                      = results.mean_In_to_Ex;
    temp=results.firings;
    meta.firings_E{iTrial}=results.firings_E;
    meta.firings_I{iTrial}=results.firings_I;
    temp2=results.time_stats.spike_rate;
    for iNeuron=1:M.p.Ne+M.p.Ni
        spike_times{iNeuron,iTrial}=temp(temp(:,2)==iNeuron,1);
        spike_rate(iNeuron,iTrial)              = temp2(iNeuron);
    end
    % Spike density functions
    SDF_temp = SDF(results.spike_raster,0.005,[1 M.p.t_max],2);
    mean_SDF_temp(:,iTrial)    = nanmean(SDF_temp.PDF_trimmed(1:M.p.Ne,:),1);
end
meta.SDF.mean_SDF_Ex            = mean_SDF_temp;
meta.SDF.grand_mean_SDF_Ex      = nanmean(mean_SDF_temp,2);
meta.SDF.grand_sem_SDF_Ex       = nansem(mean_SDF_temp,2);
meta.time_stats.spike_rate      = spike_rate';
meta.time_stats.spike_rate_mean = nanmean(meta.time_stats.spike_rate,1); meta.time_stats.spike_rate_mean(isnan(meta.time_stats.spike_rate_mean))=0.1;
meta.connectivity               = results.connectivity;
meta.spike_times                = spike_times;


meta.InExBalance.trials = meta.mean_In_to_Ex./(-1*meta.mean_Ex_to_Ex);
meta.InExBalance.trialmean = nanmean(meta.InExBalance.trials,2); 
meta.InExBalance.trialSEM  = nansem(meta.InExBalance.trials,2); 
clear iTrial results Fname temp temp2 iNeuron spike_times spike_rate
for iNeuron=1:M.p.Ne
    meta.times.count.mean(iNeuron)  = nanmean(cellfun(@numel,meta.spike_times( iNeuron,:)));
    meta.times.count.CoV(iNeuron)   = nanstd(cellfun(@numel,meta.spike_times( iNeuron,:)))./nanmean(cellfun(@numel,meta.spike_times( iNeuron,:)));
    for iTrial=1:no_trials
        temp=meta.spike_times{iNeuron,iTrial};
        try
            meta.times.first(iNeuron,iTrial)=temp(1);
        catch
            meta.times.first(iNeuron,iTrial)=NaN;
        end
    end
end

meta.times.first_mean = nanmean(meta.times.first,2);
meta.times.first_std  = nanstd(meta.times.first,2);
meta.times.first_CoV  = meta.times.first_std./meta.times.first_mean;    
meta.spectro.FFT.power_mean=nanmean(meta.spectro.FFT.power,2);meta.spectro.FFT.power_SEM=nansem(meta.spectro.FFT.power,2);
save(sprintf('%s_meta.mat',Genotype))
%% plotting
no_trials=size(meta.EIcorr.r,2);
cmap_temp=jet(no_trials);

figure; hold on
plot(meta.spectro.FFT.freq,meta.spectro.FFT.power_mean)
ciplot(meta.spectro.FFT.power_mean+meta.spectro.FFT.power_SEM,...
       meta.spectro.FFT.power_mean-meta.spectro.FFT.power_SEM,...
       meta.spectro.FFT.freq,'b')
xlabel('Frequency (Hz)'); ylabel('LFP power (dB uV^2/Hz)')
axis ([0 200 -Inf Inf])
set(gca,'XScale','log');%set(gca,'YScale','log')

figure
	set(gcf,'position',[919,5,762,215],'Units','Pixels','Color','w');
    subplot(1,2,1); hold on
    for iTrial=1:no_trials
        plot(meta.dynamics.SpikeDensity.mean_PDF_ex(1:M.p.t_max-5,iTrial),...
             meta.dynamics.SpikeDensity.mean_PDF_ex(6:M.p.t_max,iTrial),...
        'LineWidth',1.5,...
        'Color','b')
    end
    xlabel('% of Ex population (t)')
    ylabel('% of Ex population (t+5)')
    axis([0 0.3 0 0.3])
subplot(1,2,2); hold on
    for iTrial=1:no_trials
        plot(meta.dynamics.SpikeDensity.mean_PDF_ex(1:M.p.t_max,iTrial),...
             meta.dynamics.SpikeDensity.mean_PDF_in(1:M.p.t_max,iTrial),...
        'LineWidth',1.5,...
        'Color','b')%cmap_temp(iTrial,:)
    end
    xlabel('% of Ex population (t)')
    ylabel('% of In population (t)')
    axis([0 0.3 0 0.3])

figure; set(gcf,'position',[918,273,762,287],'Units','Pixels','Color','w');
subplot(1,2,1)
scatter(meta.times.first_mean(1:M.p.Ne,1)-M.input.ZAP.stim_times(1),...
        meta.times.first_std(1:M.p.Ne,1),...
        40,... %1+5*meta.time_stats.spike_rate_mean(1:M.p.Ne),... %20,... %
        meta.time_stats.spike_rate_mean(1:M.p.Ne),...
        'filled');
%         axis([0 120 0 Inf])
        xlabel('First spike latency (ms)')
        ylabel('First spike jitter (ms)')
subplot(1,2,2)
scatter(1e3*meta.connectivity.sum_inputs(1:M.p.Ne,1),...
        1e3*meta.connectivity.sum_inputs(1:M.p.Ne,2),...
            40,... %1+5*meta.time_stats.spike_rate_mean(1:M.p.Ne),... %20,... %
            meta.time_stats.spike_rate_mean(1:M.p.Ne),...
            'filled');
    axis([0 200 0 200])        
    xlabel('Net excitatory drive (nS)')
    ylabel('Net inhibitory drive (nS)')
    h=colorbar('EastOutside'); ylabel(h,'Spike count (Hz)','fontsize',12,'fontname','Helvetica');
    set(h,'Position',[0.9173,0.1568,0.0263,0.2300])
    axis square
%% Plotting - Spike density functions
figure; hold on
plot(1:M.p.t_max,meta.SDF.grand_mean_SDF_Ex,'b')
ciplot(meta.SDF.grand_mean_SDF_Ex+meta.SDF.grand_sem_SDF_Ex,...
       meta.SDF.grand_mean_SDF_Ex-meta.SDF.grand_sem_SDF_Ex,...
       1:M.p.t_max,'b')    
%% Plotting: Ex-In cross-correlation
   figure
   plot(meta.EIcorr.lags,meta.EIcorr.r)
   axis([-50 50 -Inf Inf])
%% Plotting: Spike raster
figure;
for iTrial=1:no_trials
subplot(2,5,iTrial); hold on
    plot(meta.firings_E{iTrial}(:,1),meta.firings_E{iTrial}(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',5);
    plot(meta.firings_I{iTrial}(:,1),meta.firings_I{iTrial}(:,2),'.','color',[0.5 0.8 0.9],'MarkerSize',5);
     ylabel('Neuron no.');
    axis([-90 990 0 M.p.Ne+M.p.Ni])%M.p.t_max
    axis square
end
%% Plotting: Input currents
figure
hold on
    subplot(2,1,1);hold on
    plot(M.p.tb,meta.mean_Ex_to_Ex,'LineWidth',1.2,'Color',[0.6 0.6 0.6])
    plot(M.p.tb,meta.mean_In_to_Ex,'LineWidth',1.2,'Color',[0.6 0.6 0.6])
    xlabel('Time (ms)')
    ylabel('Average input current') 
    subplot(2,1,2); hold on
    plot(M.p.tb,meta.InExBalance.trialmean,'LineWidth',1.2,'Color',[0.6 0.6 0.6])
    tempa=meta.InExBalance.trialmean;tempa(isnan(tempa))=0;
    tempb=meta.InExBalance.trialSEM;tempb(isnan(tempb))=0;
    ciplot(tempa-tempb,tempa+tempb,M.p.tb,[0.6 0.6 0.6])
    plot(M.p.tb,meta.InExBalance.trialmean,'LineWidth',1.2,'Color',[0.6 0.6 0.6])
    xlabel('Time (ms)')
    ylabel('Average In/Ex balance') 
    
