function Ssyn = DelayedSumOfExponentialsSynapse(t,t_elaps,delay,tau_rise,tau_decay)
% Simulation code to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% t         - simulation time
% t_elaps   - time since start of synaptic event
% tau_rise  - time constant for synaptic function rise 
% tau_decay - time constant for synaptic function decay
% delay     - synaptic transmission delay
%
% Outputs:
% Ssyn      - synaptic conductance calculated for each value supplied for t
% 
% Usage:
% This function approximates the conductance waveform of a synapse based on
% the sum of two exponentials. Peak synaptic amplitude is normallized to 1.
% Heavyside function negates effects of negative time. Includes terms for 
% transmission latency and elapsed time since firing.
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

% t_elaps=10
% t=-10:1:100
% 
% tau_decay=10
% tau_rise=0.5
% delay=1

t         = repmat(t,size(t_elaps));
delay     = repmat(delay,size(t_elaps));
tau_decay = repmat(tau_decay,size(t_elaps));
tau_rise  = repmat(tau_rise,size(t_elaps));

% after Roth and van Rossum
t_peak = t_elaps+delay+((tau_decay.*tau_rise)./(tau_decay-tau_rise)).*log(tau_decay./tau_rise);
f      = 1./(-exp(-(t_peak-delay-t_elaps)./tau_rise ) +  exp(-(t_peak-delay-t_elaps)./tau_decay));
Ssyn   = f.*(exp(-(delay-t_elaps)./tau_decay) - exp(-(delay-t_elaps)./tau_rise));
Ssyn(Ssyn<0)=0;



% or... after Sterrat et al.
% t_peak = t_elaps + delay + (tau_rise.*tau_decay)*log(tau_decay./tau_rise)./(tau_decay-tau_rise);
% f      = (tau_rise.*tau_decay).*(exp(-t_peak-delay/tau_decay)-exp(-t_peak./tau_rise))./(tau_decay-tau_rise);
% Ssyn   = (tau_rise.*tau_decay)*(exp(-(t-delay-t_elaps)./tau_decay)-exp(-(t-delay-t_elaps)/tau_rise))/((tau_decay-tau_rise)*f);
% Ssyn(Ssyn<0)=0;

% figure; plot(t,Ssyn)