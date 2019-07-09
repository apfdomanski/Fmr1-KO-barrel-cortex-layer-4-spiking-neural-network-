function [sig,crit] = permtest2vec(x,y,nrep,p)
% Aleks Domanski UoB 2016
%%%% permutation test for different between two samples of vectors
% Inputs:
% x and y are the N-Dim timeseries to compare [no timepoints x N repeats]
% nrep is the no. bootstraps to draw
% p is the critical pValue to report against

%Outputs:
% sig is a binary vector of same length as input, signifncalty differnt Y/N for each point in timeseries
% crit is the lower and upper bound for significant vertor separation

% Example:

% tb = 1:1000; N = 10;
% x = (randn(length(tb),N));
% y = randn(length(tb),N);
% y(200:300,:) = y(200:300,:)+2;
% 
% [p,crit] = permtest2vec(x,y,1000,0.05);
% 
% figure; hold on
% xM = mean(x,2);
% yM = mean(y,2);
% xE = std(x,[],2)./sqrt(N);
% yE = std(y,[],2)./sqrt(N);
% 
% plot(tb,xM-xE,'b',tb,xM+xE,'b')
% plot(tb,yM-yE,'r',tb,yM+yE,'r')
% 
% a = nan(size(p));a(p)=1;
% % a(x==y)=NaN;
% plot(tb,3*a,'k','LineWidth',4)
% % plot(tb,crit,':k')


if nargin<4 || isempty(p)
    p=0.05;
end
sx = size(x);
sy = size(y);

l = sx(1);
nx = sx(2);
ny = sy(2);

mnx = mean(x,2);
mny = mean(y,2);

alldata = [x y];
nall = nx+ny;

mnxsmp = zeros(nrep,l);
mnysmp = zeros(nrep,l);
% allinds = 1:nall;

ix = ceil(rand(nrep,nx)*nall);
iy = ceil(rand(nrep,ny)*nall);

for n = 1:nrep    
    cix = ix(n,:);            ciy = iy(n,:);
    cx = alldata(:,cix);      cy = alldata(:,ciy);
    mnxsmp(n,:) = mean(cx,2); mnysmp(n,:) = mean(cy,2);
end


%  "l" is the length of the vector, thus (p*100)/l is Bonferoni correction of p-value
crit = prctile(mnysmp-mnxsmp,[(p*100)/l 100-(p*100)/l],1);
sig = mny-mnx <= crit(1,:)' | mny-mnx >= crit(2,:)';
sig(mny==mnx) = 0;
end

