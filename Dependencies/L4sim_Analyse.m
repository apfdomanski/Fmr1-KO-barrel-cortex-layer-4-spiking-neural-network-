function results=L4sim_Analyse(M,results)
% Data processing for Layer 4 model results for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% results                   - single trial simulation results matrix with appended results analysis fields
%
% Usage:
% This function analyses sigle-trial simulation results and extracts statistics on Excitation/Inhibition balance in the network pool:
% - Population Ex/In balance:  Fraction of firing cells and smoothed spike density functions for each pool. 
% - Single cell Ex/In balance: Balance in input currents. 
% - Correlation between In and Ex firing pools
% - Power spectral density LFP is simulated by taking multi-taper estimates of low-pass filtered summed input currents.
% - Summary stats on spike time and rate
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

t_max    = M.p.t_max;
dt       = M.p.dt;
no_steps = t_max/dt;
Ne=M.p.Ne; Ni=M.p.Ni; N_total=Ne+Ni;
%% Excitation/Inhibition: Fraction of neurons active at each time-step
% Active populations
results.fraction_fired_E = zeros(t_max,1);
results.fraction_fired_I = zeros(t_max,1);

for t=1:t_max
    if ~isempty(results.firings_E)
    results.fraction_fired_E(t) = numel(find(results.firings_E(:,1)>=(t-1) & results.firings_E(:,1)<t));
    end
    if ~isempty(results.firings_E)
        results.fraction_fired_I(t) = numel(find(results.firings_I(:,1)>=(t-1) & results.firings_I(:,1)<t));
    end
end; clear t_idx

results.fraction_fired_E = results.fraction_fired_E./M.p.Ne;
results.fraction_fired_I = results.fraction_fired_I./M.p.Ni;
%% Spike raster and smoothed spike density functions

temp_raster= zeros(N_total,t_max);
for t=1:t_max
    temp_raster(:,t) = sum(results.spike_raster(:,M.p.tb>(t-1)&M.p.tb<=t),2)~=0;
end
    
% Spike density functions, downsampled to 1kHz
results.SDF=[];results.SDF = SDF(temp_raster,0.001,[1 M.p.t_max],1);
% results.SDF.PDF_trimmed(results.SDF.PDF_trimmed==0)=NaN;
results.SDF.mean_SDF_Ex    = nanmean(results.SDF.PDF_trimmed(1:M.p.Ne,:),1);
results.SDF.sem_SDF_Ex     = nansem(results.SDF.PDF_trimmed(1:M.p.Ne,:),1);
results.SDF.mean_SDF_In    = nanmean(results.SDF.PDF_trimmed(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
results.SDF.sem_SDF_In     = nansem(results.SDF.PDF_trimmed(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
%% Excitation/Inhibition: single cell current balance
results.mean_Ex_input_to_Ex_neurons=nanmean(results.Input_Ex(1:M.p.Ne,:),1);
results.mean_In_input_to_Ex_neurons=nanmean(results.Input_In(1:M.p.Ne,:),1);
results.mean_Ex_input_to_In_neurons=nanmean(results.Input_Ex(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
results.mean_In_input_to_In_neurons=nanmean(results.Input_In(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);

results.I_In_Ex_balance_single_cell=results.Input_In./results.Input_Ex;
results.I_In_Ex_balance_mean_Ex_neurons=nanmean(results.I_In_Ex_balance_single_cell(1:M.p.Ne,:),1);
results.I_In_Ex_balance_mean_In_neurons=nanmean(results.I_In_Ex_balance_single_cell(M.p.Ne+1:M.p.Ne+M.p.Ni,:),1);
%% Simulated LFP
% LFP is approximated as Low-pass filtered sum of all input currents (see Cavalleri...Mazoni 2014)
% Aproximate low-pass filter affect with Savitzky-Golay FIR filter
% results.LFP = 100*mean(results.I_in_total(:,:)./repmat(M.p.rm,1,M.p.t_max),1); % ex and in
results.LFP = sgolayfilt(100*nanmean(results.I_in_total(:,:)./repmat(M.p.rm,1,no_steps),1),0,5); % ex and in
% results.LFP = sgolayfilt(100*mean(results.I_in_total(1:M.p.Ne,:)./repmat(M.p.rm(1:M.p.Ne),1,M.p.t_max),1),0,5);
% results.LFP = 100*mean(results.I_in_total(1:M.p.Ne,:)./repmat(M.p.rm(1:M.p.Ne),1,M.p.t_max),1);
%% Power spectrum of simulated LFP
results.spectro = [];
fft_win         = [1 t_max-1];
fft_win= [find(M.p.tb==fft_win(1)): find(M.p.tb==fft_win(2))];
    
results.spectro.FFT.fpass      =   [1 200];
results.spectro.FFT.tapers    =   [1 1]; % time*bandwidth product and no. Slepian tapers to use
results.spectro.FFT.Fs        =   1000/dt;
results.spectro.FFT.pad       =   3;
results.spectro.FFT.trialave  =   0;


try
    %% calculate FFT power spectrum
    [results.spectro.FFT.power, results.spectro.FFT.freq] = mtspectrumc(results.LFP(fft_win),results.spectro.FFT);

    results.spectro.FFT.power = 10*log10(abs(results.spectro.FFT.power));  % Approximate power
    clear fft_win
    %% Spectrogram of simulated LFP.. see Cavalleri...Mazoni 2014
    results.spectro.analysis_window   = [1, 1 ;M.p.t_max,200]; %[min time , min freq ; max time ; max freq]
    results.spectro.p.fpass      = [0 200]; % passband in Hz
    results.spectro.p.tapers     = [1 1];
    results.spectro.p.Fs         = 1000/dt;
    results.spectro.p.pad        = 3;
    results.spectro.p.err        = [2 1];
    results.spectro.p.trialave   = 0;
    results.spectro.p.movingwin  = [0.1 0.001];

    [results.spectro.S,...
     results.spectro.T,...
     results.spectro.F,...    
     results.spectro.Serr]            = mtspecgramc(results.LFP,...
                                                     results.spectro.p.movingwin,...
                                                     results.spectro.p);
                                                 
   	 results.spectro.P                = 10*log10(abs(results.spectro.S)); % Approximate power   
end
%% Spike density estimate, downsampled to 1kHz

    results.dynamics.SpikeDensity             = SDF(temp_raster,0.001);
    % Average, SEM of Ex neurons
    results.dynamics.SpikeDensity.mean_PDF_ex = nanmean(results.dynamics.SpikeDensity.PDF_trimmed(1:M.p.Ne,:));
    results.dynamics.SpikeDensity.sem_PDF_ex  = nansem(results.dynamics.SpikeDensity.PDF_trimmed(1:M.p.Ne,:));
    % Average, SEM of In neurons
    results.dynamics.SpikeDensity.mean_PDF_in = nanmean(results.dynamics.SpikeDensity.PDF_trimmed(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
    results.dynamics.SpikeDensity.sem_PDF_in  = nansem(results.dynamics.SpikeDensity.PDF_trimmed(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
%% Net connectivity
    % N.B. connections are in format S(in,out)
    for n_ID=1:M.p.Ne+M.p.Ni
        %[Ex,In]
        results.connectivity.sum_inputs(n_ID,1)  = nansum(M.net.S(n_ID,1:M.p.Ne));
        results.connectivity.sum_inputs(n_ID,2)  = nansum(M.net.S(n_ID,M.p.Ne+1:M.p.Ne+M.p.Ni));
        temp=M.net.S(n_ID,1:M.p.Ne);
        results.connectivity.mean_inputs(n_ID,1) = nanmean(temp(temp>0));
        temp=M.net.S(n_ID,M.p.Ne+1:M.p.Ne+M.p.Ni);
        results.connectivity.mean_inputs(n_ID,2) = nanmean(temp(temp>0));
    end    
%% Spike time stats + inter-spike-interval    
    results.time_stats=[];
    results.time_stats.hist_bins=0:1:M.p.t_max;
    for n_ID=1:M.p.Ne+M.p.Ni
         results.time_stats.spike_times{n_ID} =[];
         results.time_stats.ISI{n_ID} =[];
        if ~isempty(results.firings)    
            results.time_stats.spike_times{n_ID} = results.firings(results.firings(:,2)==n_ID,1);
            results.time_stats.ISI{n_ID}         = diff(results.time_stats.spike_times{n_ID});
        end
        if ~isempty(results.time_stats.ISI{n_ID})
            results.time_stats.ISI_hist(n_ID,1:numel(results.time_stats.hist_bins))=histc(results.time_stats.ISI{n_ID},results.time_stats.hist_bins);
        else
            results.time_stats.ISI_hist(n_ID,1:numel(results.time_stats.hist_bins))=NaN;
        end
            results.time_stats.spike_rate(n_ID)=numel(results.time_stats.spike_times{n_ID})./(M.p.t_max/1000); 
    end
    results.time_stats.spike_rate_mean_Ex = nanmean(results.time_stats.ISI_hist(1:M.p.Ne,:));
    results.time_stats.spike_rate_SEM_Ex  = nansem(results.time_stats.ISI_hist(1:M.p.Ne,:));
    results.time_stats.spike_rate_mean_In = nanmean(results.time_stats.ISI_hist(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
    results.time_stats.spike_rate_SEM_In  = nansem(results.time_stats.ISI_hist(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
%% Excitation/Inhibition balance
    results.InExBalance.raw          = -1*results.Input_In./results.Input_Ex;
    results.InExBalance.mean_onto_Ex = nanmean(results.InExBalance.raw(1:M.p.Ne,:));
    results.InExBalance.SEM_onto_Ex  = nansem(results.InExBalance.raw(1:M.p.Ne,:));
    
    results.InExBalance.mean_onto_In = nanmean(results.InExBalance.raw(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
    results.InExBalance.SEM_onto_In  = nansem(results.InExBalance.raw(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
    
    results.InExBalance.SteadyState  = -1*results.connectivity.sum_inputs(:,2)./results.connectivity.sum_inputs(:,1);
    results.InExBalance.norm         = results.InExBalance.raw./repmat(results.InExBalance.SteadyState,1,no_steps);
    
    results.mean_Ex_to_Ex            = nanmean(results.Input_Ex(1:M.p.Ne,:));
    results.mean_Ex_to_In            = nanmean(results.Input_Ex(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
    
    results.mean_In_to_Ex            = nanmean(results.Input_In(1:M.p.Ne,:));
    results.mean_In_to_In            = nanmean(results.Input_In(M.p.Ne+1:M.p.Ne+M.p.Ni,:));
%% Excitation/Inhibition cross-correlation between input currents
results.EIcorr=[];

% Xcorr- population average EPSCs and IPSCs
EPSC_temp=zscore(results.mean_Ex_to_Ex); 
IPSC_temp=zscore(-1*results.mean_In_to_Ex);
[results.EIcorr.r,results.EIcorr.lags]    = xcorr(EPSC_temp(1:no_steps),IPSC_temp(1:no_steps),'coeff');
results.EIcorr.lags=results.EIcorr.lags*dt;

% Xcorr- individual cells' EPSCs and IPSCs
for n_ID=1:M.p.Ne
[results.EIcorr.r_all(:,n_ID),~]    = xcorr(zscore(results.Input_Ex(n_ID,1:no_steps)),...
                                            zscore(-1*results.Input_In(n_ID,1:no_steps)),'coeff');
end

% Xcorr - active populations
% [results.EIcorr.r,results.EIcorr.lags]    = xcorr(results.fraction_fired_E,results.fraction_fired_I,'coeff');
end
    
    