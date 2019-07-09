function ClassificationError=DecodingError(InputVector,Classes,reg)
% Leave-one out classification strength (inear discrimination error)
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Input Variables:
% InputVector:  Matrix of continuous variables [time x no. trials]
% Classes:      Event category labels (either 1 or 2) [no. trials]
% reg:          Regularization parameter
%
% Output Variables:
% ClassificationError:  Decoding strength 
%
if nargin<3, reg=0; end;
% b=waitbar(0,'Decoding...');

nTrials=length(Classes);
TrialLength=round(size(InputVector,1)/nTrials);
ClassificationError=zeros(1,TrialLength);
for t=1:TrialLength % Loop across each time bin
    %     waitbar(t/TrialLength,b,sprintf('Decoding timestep: %d / %d\n',t,TrialLength))
    X=InputVector(t:TrialLength:end,:);
    if sum(sum(X))==0
        ClassificationError(t) = nTrials;
    else
        for i=1:nTrials % Loop across each trial
            % Leave one trial out
            TrialCondition_i=Classes([1:i-1 i+1:end]);                              % Event outcome of witheld trial
            Xi=X([1:i-1 i+1:end],:);                                                % input value for witheld trial
            % Sort remaining N-1 responses into condition 1 or 2    
            X1=Xi(TrialCondition_i==1,:); k1= ~isnan(X1(:,1)); X1=X1(k1,:);  
            X2=Xi(TrialCondition_i==2,:); k2= ~isnan(X2(:,1)); X2=X2(k2,:);
            n1=length(k1); n2=length(k2);                                           % No of each type of condition
            
            PooledCovariance=((n1-1)*cov(X1)+(n2-1)*cov(X2))./(n1+n2-2);            % Pooled covariance matrix
            PooledCovariance_reg=PooledCovariance+reg*eye(size(PooledCovariance));  % Regularized pooled covariance matrix
            
            % Mahalanobis distance from each class mean at this time point
            m1=mean(X1); m2=mean(X2);                                               % Mean at each time point
            d1=(X(i,:)-m1)*PooledCovariance_reg^-1*(X(i,:)-m1)';                    % to class 1 mean
            d2=(X(i,:)-m2)*PooledCovariance_reg^-1*(X(i,:)-m2)';                    % to class 2 mean
            
            % Which is best fit?
            PredictedCondition=(d1>d2)+1;                                           % Predicted class at this time point
            
            % Accumulate classification error at this time point
            ClassificationError(t)=ClassificationError(t)+(PredictedCondition~=Classes(i));                         
        end;
    end
end;
% Normalise classification error by trial count
ClassificationError=ClassificationError./nTrials;                                                   % normalise byt no. trials


