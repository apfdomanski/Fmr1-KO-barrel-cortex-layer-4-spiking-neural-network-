function M=L4sim_MakeCurentInjectionInput(M,t_max,dt,Amp,StimDuration,plot_online)

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
% Amp                       - [1x2] [Ex pool, In pool] Amplitude of simulated current step injection (pA)
% StimDuration              - [1x2] [Ex pool, In pool] Duration of current step injection (ms)
% plot_online               - Plot the stimulation pattern results
%
% Outputs:
% M                         - Model design structure
%
% Usage:
% This script simulates a current clamp protocol in which a transient current step is
% injected into the soma of a whole-cell patched neuron.
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
M.p.tb    = (0:dt:t_max-dt); 

M.input.ZAP.A        = Amp;   %5;          %Amplitude scaling 
M.input.ZAP.Fs       = 1000/dt;            %sampling frequency (Hz)

M.input.ZAP.stim_times=0.2*t_max;

M.input.ZAP.stim_amp = [1;1];% [E,I]
                           
M.input.ZAP.stim_duration = StimDuration*M.p.dt;

%   figure; plot(x,y)

external_input_Ex   = zeros(1,length(M.p.tb));
external_input_In   = zeros(1,length(M.p.tb));

external_input_Ex(M.input.ZAP.stim_times:(M.input.ZAP.stim_times+StimDuration)/M.p.dt)=1;
external_input_In(M.input.ZAP.stim_times:(M.input.ZAP.stim_times+StimDuration)/M.p.dt)=1;


M.input.external_input(1:M.p.Ne,:)=repmat(external_input_Ex*M.input.ZAP.A(1),M.p.Ne,1);
M.input.external_input(M.p.Ne+1:M.p.Ne+M.p.Ni,:) =repmat(external_input_In*M.input.ZAP.A(2),M.p.Ni,1);
M.input.external_input=sparse(M.input.external_input);

% M.input.external_input(1:M.p.Ne,:) = M.input.external_input(1:M.p.Ne,:).*repmat(rand(M.p.Ne,1),1,length(M.p.tb));

if plot_online==1;
    figure('name','input waveform')
    plot(M.p.tb,M.input.external_input);
    ylabel('Input current (pA)') % right y-axis
    box off
end
