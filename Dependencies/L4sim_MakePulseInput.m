function M = L4sim_MakePulseInput(M,t_max,dt,Amp,StimPattern,plot_online)

% External input builder and simulation design for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% M                         - Model design structure from L4sim_DesignNetwork.m
% t_max                     - Simulation duration (ms)
% dt                        - simulation integration timestep (ms)
% Amp                       - [1x2] [Ex pool, In pool] Strength of input to Excitatory and Inhibitory neuron pools
% StimPattern               - [1x2] [No. stims, Stimulation frequency (Hz)]
% plot_online               - Plot the stimulation pattern results
%
% Outputs:
% M                         - Model design structure
%
% Usage:
% This function specifies the simulation parameters and builds the external
% pulse input structure for the thalamocortical pulse-response simulation.
% Stimulation parameters are specified as independently tunable rhythmic 
% Dirac deltas to each of the Ex and In pools.
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

if nargin<5
    plot_online=1;
end
w = warning ('off','all');

M.p.t_max = t_max;
M.p.dt    = dt;
M.p.tb    = (0:dt:t_max-dt);               % Simulation timebase (ms)

M.input.ZAP.A        = Amp;                % Amplitude scaling 
M.input.ZAP.Fs       = 1000/dt;            % sampling frequency (Hz)

firstStim_ms=10;
noStims=StimPattern(1);
ISI=round((1e3/StimPattern(2)));           %Inter-stimulus interval
for iStim=1:noStims
    M.input.ZAP.stim_times(iStim)=firstStim_ms+(iStim-1)*ISI;
end
M.input.ZAP.stim_times(end+1) = 200;

% Alternative hard coded stim times:
% M.input.ZAP.stim_times = [10 210 410 610]; % 5Hz
% M.input.ZAP.stim_times = [10 110 210 310]; % 10Hz
% M.input.ZAP.stim_times = [10 60 110 160 210]; % 20Hz
% M.input.ZAP.stim_times = [10 30 50 70 90]; % 50Hz


M.input.ZAP.stim_amp = ones(2,length(M.input.ZAP.stim_times));
M.input.ZAP.stim_duration = 2*M.p.dt;


external_input_Ex   = zeros(1,length(M.p.tb));
external_input_In   = zeros(1,length(M.p.tb));
 for stim_id=1:numel(M.input.ZAP.stim_times)
     external_input_Ex(find(M.p.tb==M.input.ZAP.stim_times(stim_id)):...
                       find(M.p.tb==M.input.ZAP.stim_times(stim_id))+M.input.ZAP.stim_duration.*M.p.dt)...
                       =M.input.ZAP.stim_amp(1,stim_id);
     external_input_In(find(M.p.tb==M.input.ZAP.stim_times(stim_id)):...
                       find(M.p.tb==M.input.ZAP.stim_times(stim_id))+M.input.ZAP.stim_duration.*M.p.dt)...
                       =M.input.ZAP.stim_amp(2,stim_id);
 end
 
% Optional: convolve input Diracs with EPSC model (exponential decay)
% ... This is redundant since the thalamocortical synaptic waveforms are
% calculated at the conductance simulation step
% x = 0:0.1:5;y = 5*x.*exp(-3*x); y=y/max(y); figure; plot(x,y)
% external_input_Ex = conv(external_input_Ex,y);
% external_input_In = conv(external_input_In,y);

external_input_Ex = external_input_Ex(1:length(M.p.tb));
external_input_Ex = external_input_Ex /max(external_input_Ex);

external_input_In = external_input_In(1:length(M.p.tb));
external_input_In = external_input_In /max(external_input_In);

M.input.external_input(1:M.p.Ne,:) = repmat(external_input_Ex*M.input.ZAP.A(1),M.p.Ne,1);
M.input.external_input(M.p.Ne+1:M.p.Ne+M.p.Ni,:) = repmat(external_input_In*M.input.ZAP.A(2),M.p.Ni,1);
M.input.external_input                           = sparse(M.input.external_input);

if plot_online==1
    figure('name','input waveform')
    plot(M.p.tb,M.input.external_input);
    ylabel('Input current (pA)') 
    box off
end
