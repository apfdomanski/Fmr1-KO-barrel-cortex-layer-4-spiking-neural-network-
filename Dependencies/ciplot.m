function fillhandle=ciplot(lower,upper,x,color,alpha,axh)
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<3 || isempty(x)
    x=1:length(lower);
end
if nargin<4 || isempty(color)
    colour='b';
end
if nargin<5 || isempty(alpha)
    alpha = 0.2;
end

if nargin<6 || isempty(axh)
   axh = gca;
end
axes(axh);
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

fillhandle= fill([x fliplr(x)],[upper fliplr(lower)],color);
    set(fillhandle,'EdgeColor',color,'FaceAlpha',alpha,'EdgeAlpha',alpha);%set edge color
end

