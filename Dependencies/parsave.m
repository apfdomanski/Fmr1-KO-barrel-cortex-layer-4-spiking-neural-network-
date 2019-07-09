function parsave(varargin)
    fname=varargin{1};    
    for i=2:nargin
       eval([inputname(i),'=varargin{i};']);  
       if i==2
            save('-mat',fname,inputname(i));
       else
           save('-mat',fname,inputname(i),'-append');
       end        
    end