function SpikeDensity = SDF(spike_raster,kernal_sigma,time_window,binwidth)
% (c) Aleksander Domanski 2009
% Convovles spike trains with Gaussian kernal
% Takes MxN spike raster, where M is cell number and N is time


%% Check input arguments
SpikeDensity=[];
if nargin<2
    kernal_sigma = 0.001;
end
if nargin<3
    time_window =[0 size(spike_raster,2)];
    % time_window =[0 1000];
end
if nargin<4
    binwidth = 1;
end


%%
% discretize time stamps to 1ms resolution.
if unique(spike_raster) == [0;1] % (i.e. pre-timestamped)
    SpikeDensity.spike_raster_stamped = spike_raster;
else
    SpikeDensity.spike_raster_stamped=zeros(size(spike_raster,1),time_window(2)-time_window(1));
    SpikeDensity.spike_raster_stamped=spike_counts(spike_raster,time_window,binwidth);
end
% convolve time stamps for each trial with Gaussian kernal of sigma=x ms

% prepare a kernel: tandard deviation = sigma ms
kernal_shoulder=3; % Time ranges form -3*st. dev. to 3*st. dev. -> make edge vector
SpikeDensity.edges= -kernal_shoulder*kernal_sigma:0.001:kernal_shoulder*kernal_sigma;
% SpikeDensity.edges= -kernal_shoulder*kernal_sigma:kernal_sigma:kernal_shoulder*kernal_sigma;

% Evaluate the Gaussian kernel
SpikeDensity.kernel =   zeros(size(SpikeDensity.edges));
SpikeDensity.kernel =   normpdf(SpikeDensity.edges, 0, kernal_sigma);

% Multiply by bin width so the probabilities sum to 1
SpikeDensity.kernel = SpikeDensity.kernel*binwidth*1E-3;


% Convolve time-stamped spike data with the kernel
SpikeDensity.PDF_trimmed=[];
for cell_id=1:size(SpikeDensity.spike_raster_stamped,1);
    % if using raw time data input
    %     SpikeDensity.PDF_trimmed(cell_id,:)=conv(SpikeDensity.spike_raster_stamped(:,cell_id),SpikeDensity.kernel);
    % if using pre-stamped data input
    SpikeDensity.PDF_trimmed(cell_id,:)=conv(spike_raster(cell_id,:),SpikeDensity.kernel);
end

% Trim out the relevant portion of the spike density result
% Find the index of the kernel center
SpikeDensity.kernal_center = ceil(length(SpikeDensity.edges)/2);

% trim outliers
SpikeDensity.PDF_trimmed=SpikeDensity.PDF_trimmed(:,...
                        SpikeDensity.kernal_center:time_window(2)+SpikeDensity.kernal_center-1);

%% clean up
SpikeDensity.kernal_sigma=kernal_sigma;
SpikeDensity.time_window=time_window;
SpikeDensity.binwidth=binwidth;
