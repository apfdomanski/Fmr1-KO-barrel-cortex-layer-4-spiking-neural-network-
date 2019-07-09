function [results] = L4sim_RunCurrentInjection(M,plot_runtime,STP_on)

% Runtime simulation for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
% This version only integrates current from a simulated whole-cell recording:
% I.E. only external current is evaluated, lacks terms for synaptic connectivity between cells.
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
% M            - Model design structure passed from L4sim_DesignNetwork.m
% plot_runtime - Options for plotting results in the loop:
%                0 = headless
%                1 = plot input vs. simulated LFP
%                2 = plot realtime synaptic weights
%                3 = plot Vm
% STP_on       - (Boolean) Include short-term plasticity calculations 
%
% Outputs:
% results      - Structure containing simulation results for further analysis
%
% Usage:
% This function runs a conductance-based spiking network simulation using predefined parameters for network connectivity and synapses. 
% Choice between leaky I&F and Izhikevich model neurons can be selected,short-term plasticity can be in/excluded and in-the-loop 
% plotting can be configured based on input arguement switches.
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
if nargin<2
    plot_runtime=0;
    % 0 = headless
    % 1 = plot input vs. simulated LFP
    % 2 = plot realtime synaptic weights
    % 3 = plot Vm
    % plot in realtime
end
switch plot_runtime
    case 1
        figure('name','simulated LFP'); hold on
    case 2
        figure('name','online weight change'); colormap((hot)); caxis([0 1])
    case 3
        figure('name','Vm');
    otherwise
end
isIandF =1;
%% Initialize
Ne=M.p.Ne; Ni=M.p.Ni; N_total=Ne+Ni;
t_max    = M.p.t_max;
dt       = M.p.dt;
no_steps = t_max/dt;
v = M.p.vr;  % Initial values of v
% v = M.p.vr + 1*randn(M.p.Ne+M.p.Ni,1);% rest potential
u = 0*v;                            % Initial values of u
results.firings        = [];	% spike timings

results.spike_raster   = zeros(N_total,no_steps);
results.Input_Ex       = zeros(N_total,no_steps);
results.Input_In       = zeros(N_total,no_steps);
results.Input_External = zeros(N_total,no_steps);

% Dynamic variables - intracortical
S_prime             = M.net.S;  % online vector tracking connection weight
D1                  = ones(size(M.net.STP.d));  % Online vector tracking depression state    
D1_at_last_spike    =  D1; %state of short-term depression at the last spike time;
last_spike_time     = -Inf(N_total,1);
last_spike_elapsed  =  Inf(N_total,1);


% Dynamic variables - thalamocortical
% D1_ext                  =  ones(N_total,1);     % Online vector tracking depression state    
D1_ext                  =  ones(size(M.net.STP.d_ext));     % Online vector tracking depression state    
D1_at_last_spike_ext    =  D1_ext;              % State of short-term depression at the last spike time;
last_spike_time_ext     = -Inf(N_total,1);
last_spike_elapsed_ext  =  Inf(N_total,1);
% results.received_input = [];                          % neurons that got an external input at this timestep
%% Runtime loop
progbar = waitbar(0,'Initializing...');
for step_no=1:no_steps

    t=M.p.tb(step_no);
    waitbar(step_no/no_steps,progbar,strcat(['Simulating ',num2str(t),'/',num2str(t_max),'ms'] ))  
    %% Deal with cells that fired/ were stimulated
    % active neurons/TC connections 
    if isIandF
        results.fired = find(v>=M.p.vt);           
        % Clamp Vm to peak
        v(results.fired)                            = M.p.vpeak(results.fired); 
    else
        results.fired = find(v>=M.p.vpeak);           
        % Reset cells that fired to V_reset 
        v(results.fired)                            = M.p.c(results.fired); 
        % update adaptation variable (u)
        u(results.fired)                            = u(results.fired)+M.p.d(results.fired); 
    end
    
    % record firing times/numbers 
    if ~isempty(results.fired)
        results.firings   =   [results.firings; t+0*results.fired, results.fired]; 
    end 
    
    results.ext_fired = M.input.external_input(:,step_no);   
    % keep synaptic depression state at the time of last firing
    D1_at_last_spike(:,results.fired)           = D1(:,results.fired);   
    D1_at_last_spike_ext(results.ext_fired==1)  = D1_ext(results.ext_fired==1);     
    %% Update spike time counters/state variables
    % Internal firing
    last_spike_time(results.fired)         = t;%t-1;
    last_spike_elapsed(results.fired)      = 0;
    last_spike_elapsed(last_spike_time~=t) =-1*(t-last_spike_time(last_spike_time~=t));
    % External input firing
    last_spike_time_ext(M.input.external_input(:,step_no)==1)         = t;
    last_spike_elapsed_ext(M.input.external_input(:,step_no)==1)      = 0;
    last_spike_elapsed_ext(M.input.external_input(:,step_no)~=1) = -1*(t-last_spike_time_ext(M.input.external_input(:,step_no)~=1));

    %% calculate synaptic conductances
    I_total =  -1*M.input.external_input(:,step_no);
    %% Update Vm
    if isIandF % leaky I&F model
        % Indices of cells in refractory period
        IDs=lt(last_spike_elapsed,-1.5);
        % cells not in refractory period integrate
        v(IDs)=v(IDs)+(-1*v(IDs) + M.p.vr(IDs) - M.p.rm(IDs).*I_total(IDs))*dt./(M.p.tm(IDs));
        % clamp spikes to peak voltage to prevent Vm explosion
        v(v>M.p.vt)=M.p.vpeak(v>M.p.vt);
        % cells in refractory period sit at reset potential
        v(~IDs) = M.p.vreset(~IDs);
    else % Izikevich model
        v=v+(M.p.k.*(v-M.p.vr).*(v-M.p.vt) -u + M.p.rm.*I_total)*dt./M.p.C;
        % clamp spikes to peak voltage to prevent Vm explosion
        v(v>M.p.vpeak)=M.p.vpeak(v>M.p.vpeak);
        u=u+M.p.a.*(M.p.b.*(v-M.p.vr)-u);
    end
    %% Run STP calculations
    if STP_on==1 % Include STP into M parameters?
        
        % connections from cells that fired depress, connections from cells that didn't fire recover exponentially
        D1(:,results.fired)       =  D1(:,results.fired).*M.net.STP.d(:,results.fired);
        D1(:,last_spike_time~=t)  =  M.net.STP.recovery(M.net.STP.tau_D1(:,last_spike_time~=t),D1(:,last_spike_time~=t),10);
        
        D1_ext(results.ext_fired==1) = D1_ext(results.ext_fired==1).*M.net.STP.d_ext(results.ext_fired==1 );
        D1_ext(results.ext_fired~=1) = M.net.STP.recovery(M.net.STP.tau_D1_ext(results.ext_fired~=1),D1_ext(results.ext_fired~=1),10);
        
        % update S
        S_prime=M.net.S.*D1;
        results.S_prime_out(:,step_no)=S_prime(1,:);

    end
    %% copy ouputs into saved vectors
    % record In/Ex for each cell at each time step
    results.D1_at_last_spike(:,step_no)  = D1_at_last_spike(1);
    results.D1_ext(:,step_no)  = D1_ext(1);
    results.I_in_total(:,step_no)     = I_total;
    results.spike_raster(unique(results.fired(:,1)),step_no)=1;    
    results.V_out(:,step_no) = v;
   
end
delete(progbar)
results.firings_E=[]; results.firings_I=[];
try
    results.firings_E=results.firings(results.firings(:,2)<Ne,:);
    results.firings_I=results.firings(results.firings(:,2)>Ne,:);
end
