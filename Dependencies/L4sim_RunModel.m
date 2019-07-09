function [results] = L4sim_RunModel(M,plot_runtime,STP_on)

% Runtime simulation for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
    plot_runtime=0; % plot in realtime?

end
switch plot_runtime
    case 1
        figure('name','Spike raster and one unit'); hold on
    case 2
        figure('name','online synaptic weight change'); colormap((hot)); caxis([0 1])
    case 3
        figure('name','Membrane potential');
    case 4
        figure('name','Synaptic conductances');
    case 5
        figure('name','Recovery from short-term depression');
    otherwise
end
isIandF =1; % 1)I&F or 0) Izhikevich model neuron
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
D1_ext                  =  ones(size(M.net.STP.d_ext));     % Online vector tracking depression state    
D1_at_last_spike_ext    =  D1_ext;              % State of short-term depression at the last spike time;
last_spike_time_ext     = -Inf(N_total,1);
last_spike_elapsed_ext  =  Inf(N_total,1);
%% Runtime loop
progbar = waitbar(0,'Initializing...');
for step_no=1:no_steps
    % Update the waitbar
    t = M.p.tb(step_no);
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

    % console outputs
    %     try
    %        [t,...
    %         M.input.external_input(1,step_no),...
    %         last_spike_time_ext(1),...
    %         last_spike_elapsed_ext(1),...
    %         D1_at_last_spike_ext(1),...
    %         D1_at_last_spike(1),...
    %         D1_ext(1),...
    %         last_spike_time(1),...
    %         last_spike_elapsed(1)]
    %        [t,M.input.external_input(N_total,t),...
    %         last_spike_time_ext(N_total),...
    %         last_spike_elapsed_ext(N_total),...
    %         D1_at_last_spike_ext(N_total),...
    %         D1_ext(N_total)]
    %     end
    %% calculate synaptic conductances
    % kinetics:   [rise(Ex cells) Decay(Ex cells)  rise(In cells) Decay(In cells)]
    % transdelay: [E2E E2I I2E I2I external]

     S_synAMPA =        [DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed(1:Ne),...
                                                         M.p.delay(1),...
                                                         M.p.k_AMPA(1),...
                                                         M.p.k_AMPA(2)) ;...
                         DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed(Ne+1:N_total),...
                                                         M.p.delay(2),...
                                                         M.p.k_AMPA(3),...
                                                         M.p.k_AMPA(4))     ];

     S_synGABA =        [DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed(1:Ne),...
                                                         M.p.delay(3),...
                                                         M.p.k_GABA(1),...
                                                         M.p.k_GABA(2)) ;...
                         DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed(Ne+1:N_total),...
                                                         M.p.delay(4),...
                                                         M.p.k_GABA(3),...
                                                         M.p.k_GABA(4))];
                                                     
     S_synNMDA =        [DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed(1:Ne),...
                                                         M.p.delay(3),...
                                                         M.p.k_NMDA(1),...
                                                         M.p.k_NMDA(2)) ;...
                                                         zeros(Ni,1)];

     S_synAMPA_ext =    [DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed_ext(1:Ne),...
                                                         M.p.delay(5),...
                                                         M.p.k_AMPA_ext(1),...
                                                         M.p.k_AMPA_ext(2)) ;...
                         DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed_ext(Ne+1:N_total),...
                                                         M.p.delay(5),...
                                                         M.p.k_AMPA_ext(3),...
                                                         M.p.k_AMPA_ext(4)) ];       
    
     S_synNMDA_ext =    [DelayedSumOfExponentialsSynapse(t,...
                                                         last_spike_elapsed_ext(1:Ne),...
                                                         M.p.delay(5),...
                                                         M.p.k_NMDA_ext(1),...
                                                         M.p.k_NMDA_ext(2)) ;...
                                                         zeros(Ni,1)        ];
    % deal with results of infinintely distant spikes (i.e. never been activated)
    S_synAMPA       (isnan(S_synAMPA))      = 0;
    S_synGABA       (isnan(S_synGABA))      = 0;
    S_synAMPA_ext   (isnan(S_synAMPA_ext))  = 0;
    S_synNMDA       (isnan(S_synNMDA))      = 0;
    S_synNMDA_ext   (isnan(S_synNMDA_ext))  = 0;

    % All synapses going INTO neuron x: M.net.S(x,:) 
    % All synapses going AWAY from neuron x: M.net.S(:,x)    
    
    I_synAMPA       = sum(D1_at_last_spike(:,1:Ne)...
                      .* M.net.S(:,1:Ne)...
                      .* repmat(S_synAMPA(1:Ne),1,N_total)',2)...
                      .* (v-M.p.E_rev_AMPA);
                  
             
    I_synGABA       = sum(D1_at_last_spike(:,Ne+1:N_total)...
                      .* M.net.S(:,Ne+1:N_total)...
                      .* repmat(S_synGABA(Ne+1:N_total),1,N_total)',2)...
                      .* (v-M.p.E_rev_GABA);                
                  
    I_synAMPA_ext   =(D1_at_last_spike_ext...
                      .* M.net.S_ext...
                      .* S_synAMPA_ext...
                      .* (v-M.p.E_rev_AMPA));
    
	% NMDA activation gate follows Jahr+Stevens 1990 and Gabbiani et al 1994
    NMDA_gate       = 1./(1+exp(-M.p.NMDA_alpha.*v).*(M.p.Mg_out/M.p.NMDA_beta));

    	
        
    I_synNMDA       = sum(D1_at_last_spike(:,1:Ne)...
                      .* M.net.S_NMDA(:,1:Ne)...
                      .* repmat(S_synNMDA(1:Ne),1,N_total)',2)...
                      .* (v-M.p.E_rev_NMDA)...
                      .* NMDA_gate;
                  
    I_synNMDA_ext   = (D1_at_last_spike_ext...
                      .* M.net.S_ext_NMDA...
                      .* S_synNMDA_ext...
                      .* (v-M.p.E_rev_AMPA))...
                      .* NMDA_gate;            
                  
    I_total = I_synAMPA_ext +  I_synGABA + I_synAMPA + I_synNMDA + I_synNMDA_ext ;%+0.1*randn(N_total,1);
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
    %% Run STP calculations on input conductances
    if STP_on==1 % Include STP into M parameters?
        
        % connections from cells that fired depress, connections from cells that didn't fire recover exponentially
        D1(:,results.fired)       =  D1(:,results.fired).*M.net.STP.d(:,results.fired);
        D1(:,last_spike_time~=t)  =  M.net.STP.recovery(M.net.STP.tau_D1(:,last_spike_time~=t),D1(:,last_spike_time~=t),10);
        
        D1_ext(results.ext_fired==1) = D1_ext(results.ext_fired==1).*M.net.STP.d_ext(results.ext_fired==1 );
        D1_ext(results.ext_fired~=1) = M.net.STP.recovery(M.net.STP.tau_D1_ext(results.ext_fired~=1),D1_ext(results.ext_fired~=1),10);
        
        % update history-dependent synaptic weight matrix (S)
        S_prime=M.net.S.*D1;
        results.S_prime_out(:,step_no)=S_prime(1,:);

    end
    %% copy ouputs into saved vectors
    % record In/Ex for each cell at each time step
    results.D1_at_last_spike(:,step_no)  = D1_at_last_spike(1);
    results.D1_ext(:,step_no)  = D1_ext(1);
    
    % Conductances:
    results.S_synAMPA_ext(:,step_no)  = S_synAMPA_ext; %Ex drive from external input (AMPA)
    results.S_synNMDA_ext(:,step_no)  = S_synNMDA_ext; %Ex drive from external input (NMDA)
    results.S_synAMPA_rec(:,step_no)  = S_synAMPA;     %Ex drive from network input (AMPA)
    results.S_synNMDA_rec(:,step_no)  = S_synNMDA;     %Ex drive from network input (NMDA)
    results.S_synGABA(:,step_no)      = S_synGABA;     %in drive from network input (GABA)
    
    %Currents:
    results.I_synAMPA_ext(:,step_no)  = I_synAMPA_ext;
    results.I_synNMDA_ext(:,step_no)  = I_synNMDA_ext;
    results.I_synAMPA_rec(:,step_no)  = I_synAMPA; 
    results.I_synNMDA_rec(:,step_no)  = I_synNMDA;
    results.I_synGABA(:,step_no)      = I_synGABA;
    results.I_in_total(:,step_no)     = I_total;
    
    results.S_Input_Ex(:,step_no)     = S_synAMPA + S_synNMDA + S_synAMPA_ext + S_synNMDA_ext;
    results.S_Input_In(:,step_no)     = S_synGABA;
    results.Input_Ex(:,step_no)       = I_synAMPA + I_synNMDA + I_synAMPA_ext + I_synNMDA_ext;
    results.Input_In(:,step_no)       = I_synGABA;

    results.spike_raster(unique(results.fired(:,1)),step_no)=1;    
    %Membrane potential 
    results.V_out(:,step_no) = v;
    %% Online plotting options
    switch plot_runtime
        case 0 % Headless
        case 1 % Plot voltage waveforms
        try
            clf
            subplot(2,1,1); hold on
            % (1) plot time stamps
            temp_E=results.firings(results.firings(:,2)<=Ne,:);
            temp_I=results.firings(results.firings(:,2)>Ne,:);
            plot(temp_E(:,1),temp_E(:,2),'.','color',[0.6 0.6 0.6],'MarkerSize',8);
            plot(temp_I(:,1),temp_I(:,2),'.','color',[0.6 0.6 0.9],'MarkerSize',8);
            area((1:step_no)*dt,10*sum(results.spike_raster(1:M.p.Ne,1:step_no)),'FaceColor',[0.6 0.6 0.6])
            area((1:step_no)*dt,10*sum(results.spike_raster(M.p.Ne+1:M.p.Ne+M.p.Ni,1:step_no)),'FaceColor',[0.6 0.6 0.9])
            
            axis([0 t -100 N_total])
            subplot(2,1,2); hold on
            plot((1:step_no)*dt,(results.V_out(Ne,1:step_no)),'color',[0.6 0.6 0.6],'LineWidth',1); hold on
            axis([0 t -100 100])
        end      
        drawnow
        case 2 % Plot short-term depression: matrix
            imagesc(vertcat(D1,repmat(D1_ext,1,100)')); drawnow
        case 3 % Plot membrane potential
            try
            plot_width=20;
            clf;
            plot((1:t)*dt,results.V_out(1:100:Ne,1:t),'color',[0.3 0.3 0.3],'LineWidth',1); hold on
            plot((1:t)*dt,150+results.V_out(Ne+1:20:N_total,1:t),'color',[0.6 0.6 0.9],'LineWidth',1);%         set(gcf,'color',[0 0 0]);set(gca,'color',[0 0 0])
            axis([0, t*dt, -100, 200]), grid on
            xlabel('time (ms)')
            ylabel('Vm (mV)')
            drawnow
        end
        case 4 % Plot synaptic conductances
            plot_width=10;
            clf;
            try
                subplot(2,2,1);hold on
                plot((t-plot_width:t)*dt,-1*results.I_synAMPA_ext(1:100:Ne,t-plot_width:t),'color',[0.3 0.3 0.3],'LineWidth',1);
                plot((t-plot_width:t)*dt,-1*results.I_synNMDA_ext(1:100:Ne,t-plot_width:t),'color',[0.6 0.3 0.3],'LineWidth',1);
                plot((t-plot_width:t)*dt,   results.I_synGABA    (1:100:Ne,t-plot_width:t),'color',[0.3 0.6 0.3],'LineWidth',1);
                axis([(t-plot_width)*dt, t*dt, -1, 1]), grid on
                title('Ex cells: Input')
                xlabel('time (ms)')
                ylabel('Input Conductance (nC)')
                subplot(2,2,2); hold on
                plot((t-plot_width:t)*dt,-1*results.I_synAMPA_ext(Ne+1:20:N_total,t-plot_width:t),'color',[0.6 0.6 0.9],'LineWidth',1);
                plot((t-plot_width:t)*dt,   results.I_synGABA    (Ne+1:20:N_total,t-plot_width:t),'color',[0.3 0.6 0.3],'LineWidth',1);
                axis([(t-plot_width)*dt, t*dt, -1, 1]), grid on
                title('In cells: Input')
                xlabel('time (ms)')
                ylabel('Input Conductance (nC)')
                subplot(2,2,3); hold on
                plot((t-plot_width:t)*dt,-1*results.I_synAMPA_rec(1:100:Ne,t-plot_width:t),'color',[0.3 0.3 0.3],'LineWidth',1);
                plot((t-plot_width:t)*dt,-1*results.I_synNMDA_rec(1:100:Ne,t-plot_width:t),'color',[0.6 0.3 0.3],'LineWidth',1);
                plot((t-plot_width:t)*dt,   results.I_synGABA    (1:100:Ne,t-plot_width:t),'color',[0.3 0.6 0.3],'LineWidth',1);
                axis([(t-plot_width)*dt, t*dt, -1, 1]), grid on
                title('Ex cells: recurrent')
                xlabel('time (ms)')
                ylabel('Input Conductance (nC)')
                subplot(2,2,4)
                plot((t-plot_width:t)*dt,-1*results.I_synAMPA_rec(Ne+1:20:N_total,t-plot_width:t),'color',[0.6 0.6 0.9],'LineWidth',1);
                plot((t-plot_width:t)*dt,   results.I_synGABA    (Ne+1:20:N_total,t-plot_width:t),'color',[0.3 0.6 0.3],'LineWidth',1);
                axis([(t-plot_width)*dt, t*dt, -1, 1]), grid on
                title('In cells: recurrent')
                xlabel('time (ms)')
                ylabel('Input Conductance (nC)')
                drawnow
            end
        case 5 % Plot short-term depression: recovery
            plot_width=10;
            clf;
            try
                subplot(2,1,1);hold on; cla
                plot((1:step_no)*dt,results.D1_at_last_spike(1:step_no),'color',[0.6 0.6 0.9],'LineWidth',1);
                subplot(2,1,2);hold on;cla
                plot((1:step_no)*dt,results.D1_ext(1:step_no),'color',[0.6 0.6 0.9],'LineWidth',1);
            catch
            end
    end
end
%% Clear up
delete(progbar)
results.firings_E=[]; results.firings_I=[];
try
    results.firings_E=results.firings(results.firings(:,2)<Ne,:);
    results.firings_I=results.firings(results.firings(:,2)>Ne,:);
end
