function M=L4sim_DesignNetwork( Ne,Ni,...     
                                d_ext2I,tau_D1_ext2I,...
                                d_ext2E,tau_D1_ext2E,...
                                d_I2E  ,tau_D1_I2E,...
                                d_E2I  ,tau_D1_E2I,...
                                d_I2I  ,tau_D1_I2I,...
                                d_E2E  ,tau_D1_E2E,...
                                E2E,...
                                E2I,...
                                I2E,...
                                I2I,...
                                ext_input,...
                                plot_online,...
                                Condition)

% Network builder for L4 thalamocortical integration spiking model to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% Ne,Ni                     - number of Ex, In neurons
% d_ext2I,tau_D1_ext2I      - TC->In STP parameters: dep. per stim and recovery time constant
% d_ext2E,tau_D1_ext2E      - TC->Ex STP parameters: dep. per stim and recovery time constant
% d_I2E  ,tau_D1_I2E        - In->Ex STP parameters: dep. per stim and recovery time constant
% d_E2I  ,tau_D1_E2I        - Ex->In STP parameters: dep. per stim and recovery time constant
% d_I2I  ,tau_D1_I2I        - In->In STP parameters: dep. per stim and recovery time constant
% d_E2E  ,tau_D1_E2E        - Ex->Ex STP parameters: dep. per stim and recovery time constant
% E2E                       - Ex to Ex connection [p(con), mean strength]...
% E2I                       - Ex to In connection [p(con), mean strength]...
% I2E                       - In to Ex connection [p(con), mean strength]...
% I2I                       - In to In connection [p(con), mean strength]...
% ext_input                 - External gmax: AMPA[Ex, In], NMDA[Ex] (nS)
% plot_online               - (Boolean) Plot matrix of connectivity and summary stats
% Condition)                - (String: 'WT' or 'KO') Preset physiology values each genotype
%
% Outputs:
% M                         - Model design structure
%
% Usage:
% This function builds intrinsic parameters and a synaptic connectivity
% matrix for a recurrent spiking neural network with external stimulation
% and synapse-specific short-term plasticity.
% All synapses going INTO neuron x: M.net.S(x,:) 
% All synapses going FROM  neuron x: M.net.S(:,x)
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

%% Cell parameters and distributions

% M.p.rm        Input resistance (MO) 
% M.p.C         Membrane capacitance (pF)
% M.p.vr        Resting potential (mV)
% M.p.vreset    IandF reset potential (mV)
% M.p.vt        Firing threshold (mV)
% M.p.vpeak     Peak potential (mV)

% no. neurons                  
M.p.Ne = Ne; 
M.p.Ni = Ni; 
% random seeds 
M.p.re             = randn(M.p.Ne,1);       % For Ex neurons   	
M.p.ri             = randn(M.p.Ni,1);       % For In neurons
% optional: parameters for Izhikevich implementation
M.p.a              = [0.05*ones(M.p.Ne,1);... 
                      0.02*ones(M.p.Ni,1)];
M.p.b              = [0.2*ones(M.p.Ne,1);...
                      0.2*ones(M.p.Ni,1)];
M.p.c              = [-55*ones(M.p.Ne,1);...
                      -65*ones(M.p.Ni,1)];
M.p.d              = [10*ones(M.p.Ne,1);...
                       8*ones(M.p.Ni,1)];
M.p.k              = [0.25*ones(M.p.Ne,1);...
                       1*ones(M.p.Ni,1)];

M.p.vr             = [-65*ones(M.p.Ne,1);...  % Resting potential
                      -65*ones(M.p.Ni,1)];
    
M.p.vreset         = [-65*ones(M.p.Ne,1);...  % Spike reset potential
                      -85*ones(M.p.Ni,1)];

M.p.vt             = [-45*ones(M.p.Ne,1);...  % Threshold potential
                      -45*ones(M.p.Ni,1)];


% Genotype-specific intrinsic parameters: Cm and Rm covary with shared and private
% contributions to their variance
if strcmp(Condition, 'WT')
    M.p.rm              = [410+50*M.p.re;...
                           100+20*M.p.ri];
    M.p.C              = [90+10*M.p.re+10*randn(M.p.Ne,1);...
                          40+5*M.p.ri+10*randn(M.p.Ni,1)];                           
elseif strcmp(Condition, 'KO')
    M.p.rm             = [600+100*M.p.re;...
                          300+80*M.p.ri];
    M.p.C              = [90+10*M.p.re+15*randn(M.p.Ne,1);...
                          40+5*M.p.ri+10*randn(M.p.Ni,1)];    
end

M.p.tm             =  (M.p.rm*1e6).*(M.p.C*1e-9); % membrane tau (ms)
                           
M.p.vpeak          = [35*ones(M.p.Ne,1);... %Spike peak voltage
                      20*ones(M.p.Ni,1)];

% Jitter the rest potential
M.p.vr = M.p.vr + 1*randn(M.p.Ne+M.p.Ni,1);
%% net connection parameters
M.net.E2E      = E2E; % Ex to Ex
M.net.E2I      = E2I; % Ex to In 
M.net.I2E      = I2E; % In to Ex
M.net.I2I      = I2I; % In to In
%% Conducance-based Network parameters
%%%% Reversal potentials
M.p.E_rev_AMPA = 0;
M.p.E_rev_NMDA = 0;
M.p.E_rev_GABA = -71;
%%%% NMDA gating variables: Jahr & Stevens 1990
M.p.NMDA_alpha = 0.062; %mV^-1
M.p.NMDA_beta  = 3.57;  %mM
M.p.Mg_out     = 1.2;   %mM
%%%% Synaptic delay: [E2E, E2I, I2E, I2I, external Ex input]
M.p.delay = [0.5, 0.5, 0.5, 0.5, 0];
%%%% Synaptic kinetics: [rise(Ex cells), Decay(Ex cells),  rise(In cells), Decay(In cells)]
if strcmp(Condition, 'WT')
    M.p.k_AMPA     = [0.6 3 0.1 1];
    M.p.k_AMPA_ext = [0.6 3 0.1 1];
    M.p.k_NMDA     = [3 75];
    M.p.k_NMDA_ext = [3 75];
    M.p.k_GABA     = [0.7 8 0.7 8];
elseif strcmp(Condition, 'KO')
    M.p.k_AMPA     = [0.6 3 0.1 1];
    M.p.k_AMPA_ext = [0.6 3 0.1 1];
    M.p.k_NMDA     = [3 75];
    M.p.k_NMDA_ext = [3 75];
    M.p.k_GABA     = [1 12 1 12];
end
%% Synaptic network connection parameters
% (1) Rectified Gaussian distributions
M.p.ree             = abs(randn(M.p.Ne,M.p.Ne));          	
M.p.rei             = abs(randn(M.p.Ne,M.p.Ni));
M.p.rie             = rand(M.p.Ni,M.p.Ne);  
M.p.rii             = rand(M.p.Ni,M.p.Ni);

% (2) All to all
% M.p.ree             = ones(M.p.Ne,M.p.Ne);          	
% M.p.rie             = ones(M.p.Ni,M.p.Ne);  
% M.p.rei             = ones(M.p.Ne,M.p.Ni);
% M.p.rii             = ones(M.p.Ni,M.p.Ni);

% (3) Lognormal strength
% M.p.ree             = lognrnd(1,0.001,M.p.Ne,M.p.Ne);          	
% M.p.rei             = lognrnd(1,0.5,M.p.Ne,M.p.Ni);

%
% ./S: AMPA connections
% ./S_NMDA: NMDA connections (Only for Ex to Ex)
M.net.S       =[M.net.E2E(2)*M.p.ree , M.net.I2E(2)*M.p.rei;...
                M.net.E2I(2)*M.p.rie , M.net.I2I(2)*M.p.rii];                    
M.net.S_NMDA  =[M.net.E2E(3)*M.p.ree , zeros(M.p.Ne,M.p.Ni);...
                zeros(M.p.Ni,M.p.Ne),  zeros(M.p.Ni,M.p.Ni)];


%% Set network connection probability 
prune_E2E = ones(M.p.Ne,M.p.Ne); prune_E2E(rand(size(prune_E2E))>M.net.E2E(1))=0;
prune_I2E = ones(M.p.Ne,M.p.Ni); prune_I2E(rand(size(prune_I2E))>M.net.I2E(1))=0;
prune_E2I = ones(M.p.Ni,M.p.Ne); prune_E2I(rand(size(prune_E2I))>M.net.E2I(1))=0;
prune_I2I = ones(M.p.Ni,M.p.Ni); prune_I2I(rand(size(prune_I2I))>M.net.I2I(1))=0;

prune_matrix   = [prune_E2E,prune_I2E ; prune_E2I,prune_I2I];
M.net.S(prune_matrix~=1)=0; 
M.net.S_NMDA(prune_matrix~=1)=0; 
clear prune_E2E prune_I2E prune_E2I prune_I2I prune_matrix

% (3) remove autapses
M.net.S(logical(eye(size(M.net.S)))) = 0;
M.net.S_NMDA(logical(eye(size(M.net.S_NMDA)))) = 0;
%% External synaptic drive drive strength
% Add some jitter to the strength of the External input to Ex cells
re = 0.1*rand(M.p.Ne,1); 
ri = 0.05*rand(M.p.Ni,1);
% (1) all share same strength
M.net.S_ext = [ext_input(1)*ones(M.p.Ne,1); ext_input(2)*ones(M.p.Ni,1)];
M.net.S_ext_NMDA = [ext_input(3)*ones(M.p.Ne,1); zeros(M.p.Ni,1)];        

% (2) Subtract a random fraction to jitter
M.net.S_ext =      [M.net.S_ext(1:M.p.Ne) - re*ext_input(1);...
                    M.net.S_ext(M.p.Ne+1:M.p.Ni+M.p.Ne) - ri*ext_input(2)];
                    
M.net.S_ext_NMDA = [M.net.S_ext_NMDA(1:M.p.Ne)-re*ext_input(3);...
                    zeros(M.p.Ni,1)];      
                
% (3) Silence a random fraction of Ex input
r=[rand(M.p.Ne,1);ones(M.p.Ni,1)];
M.net.S_ext(r<0.2)=0;
M.net.S_ext_NMDA(r<0.2)=0;
%% Short-term depression parameters
% Varela... Nelson (1997) inspired STP model, single exponential recovery after constant fractional depression.    
M.net.STP.recovery   = @(tau,D,t) 1-(1-D).* exp((1-t)./tau);

% (1) jittered 
M.net.STP.d         = [d_E2E*1+0.01*randn(M.p.Ne,M.p.Ne) , d_I2E*1+0.01*randn(M.p.Ne,M.p.Ni);...
                       d_E2I*1+0.01*randn(M.p.Ni,M.p.Ne) , d_I2I*1+0.01*randn(M.p.Ni,M.p.Ni)];
M.net.STP.tau_D1    = [tau_D1_E2E*1+0.01*randn(M.p.Ne,M.p.Ne) , tau_D1_I2E*1+0.01*randn(M.p.Ne,M.p.Ni);...
                       tau_D1_E2I*1+0.01*randn(M.p.Ni,M.p.Ne) , tau_D1_I2I*1+0.01*randn(M.p.Ni,M.p.Ni)];                   
M.net.STP.d_ext      = [d_ext2E*1+0.01*randn(M.p.Ne,1); d_ext2I*1+0.01*randn(M.p.Ni,1)];
M.net.STP.tau_D1_ext = [tau_D1_ext2E*1+0.01*randn(M.p.Ne,1); tau_D1_ext2I*1+0.01*randn(M.p.Ni,1)];

% (2) all the same
% M.net.STP.d         = [d_E2E*ones(M.p.Ne,M.p.Ne) , d_I2E*ones(M.p.Ne,M.p.Ni);...
%                                d_E2I*ones(M.p.Ni,M.p.Ne) , d_I2I*ones(M.p.Ni,M.p.Ni)];
% M.net.STP.tau_D1    = [tau_D1_E2E*ones(M.p.Ne,M.p.Ne) , tau_D1_I2E*ones(M.p.Ne,M.p.Ni);...
%                                tau_D1_E2I*ones(M.p.Ni,M.p.Ne) , tau_D1_I2I*ones(M.p.Ni,M.p.Ni)];                   
% M.net.STP.d_ext      = [d_ext2E*ones(M.p.Ne,1); d_ext2I*ones(M.p.Ni,1)];
% M.net.STP.tau_D1_ext = [tau_D1_ext2E*ones(M.p.Ne,1); tau_D1_ext2I*ones(M.p.Ni,1)];
%% Calculate connection statistics
M.net.stats.S_ex_bins = (0:0.2:5)*1e-4;
M.net.stats.S_ex = histc(M.net.S(1:M.p.Ne,1:M.p.Ne),M.net.stats.S_ex_bins);               
M.net.stats.S_ex       = mean(M.net.stats.S_ex,2); 
M.net.stats.S_ex       = M.net.stats.S_ex./sum(M.net.stats.S_ex);

M.net.stats.S_in_bins = (0:0.2:5)*1e-4;
M.net.stats.S_in = histc(M.net.S(1:M.p.Ne,...
                                 M.p.Ne+1:M.p.Ne+M.p.Ni),...
                         M.net.stats.S_ex_bins);   
 M.net.stats.S_in       = mean(M.net.stats.S_in,2); 
 M.net.stats.S_in       = M.net.stats.S_in./sum(M.net.stats.S_in);

 M.net.Ex_strNames      = num2cell(1:M.p.Ne);
 M.net.In_strNames      = num2cell(M.p.Ne+1:M.p.Ne+M.p.Ni);
 M.net.All_strNames     = num2cell(1:M.p.Ne+M.p.Ni);
%% plotting options
if plot_online==1
    %% Distribution of intrinsic parameters
    figure;
    cmap = jet;
    subplot(4,1,1);hold on
    bins=0:10:1000;
    bar(bins,histc(M.p.rm(1:M.p.Ne),bins),'Facecolor',[0.6 0.6 0.6],'Edgecolor',[0.6 0.6 0.6 ]);
    bar(bins,histc(M.p.rm(M.p.Ne+1:M.p.Ne+M.p.Ni),bins),'Facecolor',[0.6 0.6 1],'Edgecolor',[0.6 0.6 1]);
    axis ([min(bins) max(bins) 0 Inf])
    xlabel('Input resistace (MOhms)')
    ylabel('No. Neurons')
    subplot(4,1,2);hold on
    bins=0:1:200;
    bar(bins,histc(M.p.C(1:M.p.Ne),bins),'Facecolor',[0.6 0.6 0.6],'Edgecolor',[0.6 0.6 0.6 ]);
    bar(bins,histc(M.p.C(M.p.Ne+1:M.p.Ne+M.p.Ni),bins),'Facecolor',[0.6 0.6 1],'Edgecolor',[0.6 0.6 1]);
    axis ([min(bins) max(bins) 0 Inf])
    xlabel('Time constant (ms)')
    ylabel('No. Neurons')
    subplot(2,1,2);hold on
    scatter(M.p.rm(1:M.p.Ne),...
        M.p.C(1:M.p.Ne),'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6 ]);
    scatter(M.p.rm(M.p.Ne+1:M.p.Ne+M.p.Ni),...
        M.p.C(M.p.Ne+1:M.p.Ne+M.p.Ni),'MarkerEdgeColor',[0.6 0.6 1],'MarkerFaceColor',[0.6 0.6 1]);
    axis ([0 1000 0 200])
    xlabel('Input resistace (M\Omega)')
    ylabel('Time constant (ms)')
    %% Distribution of intrinsic parameters
    figure('Color','w'); hold on
    h=scatterhist(M.p.rm,M.p.tm,'group',[repmat({'Ex'},M.p.Ne,1);repmat({'In'},M.p.Ni,1)],...
        'Location','Northwest','Direction','out','Color',[0.6 0.6 0.6;0.6 0.6 1],...
        'Style','bar','LineWidth',[2,2],...
        'kernel','overlay','Legend','off');
    hp=get(h(1),'children'); % handle for plot inside scaterplot axes
    set(hp,'Marker','.','MarkerSize',20);
    xlabel('Input resistace (M\Omega)')
    ylabel('Time constant (ms)')
    axis([0 1000 0 150])
    legend('Ex neurons','In neurons','Location','southoutside','Orientation','horizontal'); legend('boxoff')    
    %% Synaptic weight matrix
    figure('color','w'); hold on
    title('Synaptic weight matrix')
    hold on
    imagesc((M.net.S*1000)');
    rectangle('Position',[1 1 M.p.Ne M.p.Ne],'EdgeColor',[0.6 0.6 0.6],'LineWidth',2)
    text(M.p.Ne/2,M.p.Ne/2,'Ex to Ex','color','w','HorizontalAlignment','center','Fontsize',18)
    rectangle('Position',[M.p.Ne+1 1 M.p.Ni M.p.Ne],'EdgeColor',[0.6 0.6 0.6],'LineWidth',2)
    text(M.p.Ne/2,M.p.Ne+M.p.Ni/2,'In to Ex','color','w','HorizontalAlignment','center','Fontsize',18)
    rectangle('Position',[M.p.Ne+1 M.p.Ne+1 M.p.Ni M.p.Ni],'EdgeColor',[0.6 0.6 0.6],'LineWidth',2)
    text(M.p.Ne+M.p.Ni/2,M.p.Ne +M.p.Ni/2,'In to In','color','w','HorizontalAlignment','center','Fontsize',18)
    rectangle('Position',[1 M.p.Ne+1 M.p.Ne M.p.Ni],'EdgeColor',[0.6 0.6 0.6],'LineWidth',2)
    text(M.p.Ne+M.p.Ni/2,M.p.Ne/2,'Ex to In','color','w','HorizontalAlignment','center','Fontsize',18)
    axis square
    set(gca,'xtick',[1 M.p.Ne M.p.Ne+M.p.Ni],'ytick',[1 M.p.Ne M.p.Ne+M.p.Ni],...
        'xticklabel',[1 M.p.Ne M.p.Ni],'yticklabel',[1 M.p.Ne M.p.Ni])
    colormap((jet)); caxis([0 0.5])
    
    cb = colorbar;
    cb.Label.String = 'Synaptic Strength (|nS|)';
    cb.Label.FontSize = 12;
    xlabel('From neuron no.')
    ylabel('To neuron no.')
    axis([1 M.p.Ne+M.p.Ni+1 1 M.p.Ne+M.p.Ni+1])
    %% Distribution of synaptic strength
    figure;hold on; box off
    plot(M.net.stats.S_ex_bins*1000, mean(M.net.stats.S_ex,2),'color',[0.6 0.6 0.6],'LineWidth',2)
    xlabel('Ex to Ex input strength (nS)'); ylabel('Fraction of connections');
    set(gca,'yscale','log')
    plot(M.net.stats.S_in_bins*1000,mean(M.net.stats.S_in,2),'color',[0.5 0.8 0.9],'LineWidth',2)
    xlabel('Connection strength (|nS|)'); ylabel('Fraction of connections');
    set(gca,'yscale','log')
    axis([0 0.5 1e-6 1])
    legend('Ex synapses','In Synapses');legend('boxoff')
    %% Plot external AMPA and NMDA drive to each neuron
    figure; hold on;
    plot(M.net.S_ext,'b')
    plot(M.net.S_ext_NMDA,'g')
    xlabel('Neuron no.'); ylabel('TC input g_{max} (nS)');
    legend('AMPA','NMDA'); legend('boxoff')
end