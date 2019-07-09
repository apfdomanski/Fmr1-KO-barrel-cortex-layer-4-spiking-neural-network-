function save_parfor(fname, data)
% Asynchronous saving from within a parfor loop. Pass transparent variable names into function namespace.
% APFD UoB 2015 aleks.domanski@bristol.ac.uk
% @TODO:02/02/2017: Fix r2016b issue
% 
varName=genvarname(inputname(2));
% varName = matlab.lang.makeValidName(inputname(2));
eval([varName '=data;']);
try
    save(fname,varName,'-append');
catch
    save(fname,varName);
end
end


% varName=genvarname(inputname(2));
% eval([varName '=data;']);
% try
%     save(fname,varName,'-append');
% catch
%     save(fname,varName);
% end
% end

% function save_parfor(varargin)
% 
%    fname=varargin{1};    
%     for i=2:nargin
%        eval([inputname(i),'=varargin{i};']);  
%        if i==2
%             save('-mat',fname,inputname(i));
%        else
%            save('-mat',fname,inputname(i),'-append');
%        end        
%     end