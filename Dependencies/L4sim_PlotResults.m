function L4sim_PlotResults(M,results)
% Data plotting for Layer 4 model results for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% Inputs:
% M                         - Model design structure from L4sim_DesignNetwork and L4sim_MakePulseInput
% results                   - single trial simulation results matrix
%
% Outputs:
% N/A
%
% Usage:
% This function plots results of a single model trial:
% - Spike rasters and smoothed spike density plots
% - Simulated membrane potential and power spectra/spectrograms
% - Synaptic currents 
% - Measures of Excitation/Inhibition
% - Single cell firing summary stats
%
% N.B. frequency-domain analysis requires the Chronux toolbox, redistributed under the GNU  General Public License ver.2
% http://chronux.org/forum/ and "Observed Brain Dynamics", Partha Mitra and Hemant Bokil, Oxford University Press, New York, 2008.
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
%% Raw Vm data 
figure
subplot (2,1,1); hold on
for N_id=1:M.p.Ne
    patchline((M.p.tb)',results.V_out(N_id,:)','edgecolor',[0.6 0.6 0.6],'linewidth',1,'edgealpha',1);
end
patchline((M.p.tb)',results.V_out(1,:)','edgecolor',[0 0 0],'linewidth',2,'edgealpha',1);
axis off

subplot (2,1,2); hold on
for N_id=M.p.Ne+1:M.p.Ne+M.p.Ni
    patchline((M.p.tb)',results.V_out(N_id,:)','edgecolor',[0.5 0.8 0.9],'linewidth',1,'edgealpha',1);
end
patchline((M.p.tb)',results.V_out(M.p.Ne+1,:)','edgecolor',[0.3 0.5 0.9],'linewidth',2,'edgealpha',1);
axis off
%% Synaptic currents
tmax=M.p.t_max;
figure
set(gcf,'position',[1,510,969,500],'Units','Pixels','Color','w');
% Input waveform
subplot(3,2,1); hold on
    plot(M.p.tb,(results.S_synAMPA_ext(1,:)),'-k');
    plot(M.p.tb,(results.S_synNMDA_ext(1,:)),'-r');
    plot(M.p.tb,(results.S_synGABA(1,:)),'-g');    
    axis([0 tmax 0 1])
    legend('TC AMPA','TC NMDA','GABA')
    title('conductances')
subplot(3,2,2); hold on
    plot(M.p.tb,(results.I_synAMPA_ext(1,:)),'-k');
    plot(M.p.tb,(results.I_synNMDA_ext(1,:)),'-r');
    plot(M.p.tb,(results.I_synGABA(1,:)),'-g');    
    axis([0 tmax -1 1])
    legend('TC AMPA','TC NMDA','GABA')
    title('currents')
subplot(3,2,3); hold on
    plot(M.p.tb,(results.S_synAMPA_rec(1,:)),'-k');
    plot(M.p.tb,(results.S_synNMDA_rec(1,:)),'-r');
    axis([0 tmax 0 1])
    legend('recurrent AMPA','reccurent NMDA','GABA')
    title('conductances')
    
subplot(3,2,4); hold on
    plot(M.p.tb,(results.I_synAMPA_rec(1,:)),'-k');
    plot(M.p.tb,(results.I_synNMDA_rec(1,:)),'-r');
    plot(M.p.tb,(results.I_synGABA(1,:)),'-g');    
    axis([0 tmax -1 1])
    legend('recurrent AMPA','reccurent NMDA','GABA')
    title('currents')
%% Main figure
figure
set(gcf,'position',[1,510,969,500],'Units','Pixels','Color','w');

% spike raster + input waveform
subplot(3,1,1); hold on
    plot(results.firings_E(:,1),100+results.firings_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',5);
    plot(results.firings_I(:,1),100+results.firings_I(:,2),'.','color',[0.5 0.8 0.9],'MarkerSize',5);
    ylabel('Neuron no.');
    plot(M.p.tb,100.*(results.S_synAMPA_ext(1,:)),'-g');
    axis([0 M.p.t_max 0 100+M.p.Ne+M.p.Ni])

% fraction of co-active cells
subplot(3,1,2); hold on
    hold on
    tb=1:M.p.t_max;    
    area(tb,results.SDF.mean_SDF_In,'FaceColor',[0.5 0.8 0.9],'EdgeColor',[0.5 0.8 0.9])
    area(tb,results.SDF.mean_SDF_Ex,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6])
    axis([0 M.p.t_max 0 Inf])
    ylabel('Fraction co-active');

% plot dynamic In/Ex balance       
subplot(3,1,3); hold on
    plot(M.p.tb,-1*results.mean_Ex_to_Ex,'LineWidth',2,'Color',[0.6 0.6 0.6])
    plot(M.p.tb,-1*results.mean_In_to_Ex,'LineWidth',2,'Color',[0.6 0.6 0.6])
    xlabel('Time (ms)')
    ylabel('Average input current')        
%% Ex-In correlation
figure; hold on
set(gcf,'position',[1018,589,323,271],'Units','Pixels','Color','w'); hold on
try
    ciplot(results.EIcorr.r+results.EIcorr.r_SEM,...
    results.EIcorr.r-results.EIcorr.r_SEM,...
    results.EIcorr.lags','k')
end
plot(results.EIcorr.lags,results.EIcorr.r,'LineWidth',2,'Color','k')
line([0 0],[0 max(results.EIcorr.r)],'color','k','LineStyle',':')
axis ([-50 50 -Inf Inf])
xlabel('Time lag'); ylabel('E-I correlation'); 
%% Simulated LFP
figure; 
subplot(2,1,1);hold on

    plot(results.LFP,'k','LineWidth',2)
    xlabel('Time (ms)'); ylabel('Simulated LFP (mV)');

    axis off

try
%% LFP spectrogram
    subplot(2,1,2);hold on
    imagesc(1000*results.spectro.T,fliplr(results.spectro.F),rot90(results.spectro.P));
    axis tight;  axis xy
    xlabel('Time (ms)'); ylabel('Frequency (Hz)');
    axis([0 M.p.t_max 0 200])
    axis([-Inf Inf 0 200])
    caxis([-50 -25])

%%  LFP Power spectrum 
    figure; hold on
    plot(results.spectro.FFT.freq,smooth_hist(results.spectro.FFT.power),'color',[0.6 0.6 0.6],'LineWidth',1.5)
    xlabel('Frequency (Hz)'); ylabel('(dB uV^2/Hz)')
    axis ([0 200 -Inf Inf])
    set(gca,'XScale','log');%set(gca,'YScale','log')
end
%% Ex vs. In fraction active
figure; hold on
    plot(results.dynamics.SpikeDensity.mean_PDF_ex(1:M.p.t_max),...
         results.dynamics.SpikeDensity.mean_PDF_in(1:M.p.t_max),'color','k','LineWidth',2)
    axis_max=1.2*max(horzcat(results.dynamics.SpikeDensity.mean_PDF_ex(1:M.p.t_max),results.dynamics.SpikeDensity.mean_PDF_in(1:M.p.t_max)));
    plot([0 axis_max],[0 axis_max],':k')
    set(gcf,'position',[588,612,199,200],'Units','Pixels','Color','w');
    axis ([0 axis_max 0 axis_max]); axis square; box on       
    xlabel('Fraction Ex cells active');ylabel('Fraction In cells active');
%% ISI distribution    
    figure; hold on
    set(gcf,'position',[588,273,330,288],'Units','Pixels','Color','w');
    plot(results.time_stats.hist_bins,results.time_stats.spike_rate_mean_In,'Linewidth',2,'color',[0.5 0.8 0.9])
    plot(results.time_stats.hist_bins,results.time_stats.spike_rate_mean_Ex,'Linewidth',2,'color',[0.6 0.6 0.6])
    xlabel('ISI (ms)')
    ylabel('Count')
    set(gca,'xscale','log')
%% Relationship between firing rate and steady-state input strength
    figure; 
	set(gcf,'position',[918,273,762,287],'Units','Pixels','Color','w');

    subplot(1,2,1);hold on % mean firing rate vs steady-state Ex inputs   
    scatter(1e3*results.connectivity.sum_inputs(1:M.p.Ne,1),...
            results.time_stats.spike_rate(1:M.p.Ne),'.k')  
	scatter(1e3*results.connectivity.sum_inputs(M.p.Ne+1:M.p.Ne+M.p.Ni,1),...
            results.time_stats.spike_rate(M.p.Ne+1:M.p.Ne+M.p.Ni),'.b')  
    xlabel('Net excitatory drive (nS)')
    ylabel('Spike count (Hz)')
    legend('Ex neurons','In neurons','Location','NorthWest')
    axis([0 Inf 0 Inf])
    axis square	
    subplot(1,2,2);hold on % mean firing rate vs steady-state Ex and In inputs
    scatter(1e3*results.connectivity.sum_inputs(1:M.p.Ne,1),...
            1e3*results.connectivity.sum_inputs(1:M.p.Ne,2),...
            20,... %results.time_stats.spike_rate(1:M.p.Ne)
            results.time_stats.spike_rate(1:M.p.Ne),...
            'filled');
    xlabel('Net excitatory drive (nS)')
    ylabel('Net inhibitory drive (nS)')
    h=colorbar('EastOutside'); ylabel(h,'Spike count (Hz)','fontsize',12,'fontname','Helvetica');
    set(h,'Position',[0.9173,0.1568,0.0263,0.2300])
    axis square
end