% Simulation analysis code to accompany: Domanski,Booker,Wyllie,Isaac,Kind (2018)
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
% Usage:
% This script recapitulates the simulation results for Fig 10D-E in the above paper.
% - Spike trains from distractor stimulation trials are compared againt those
%   from the reference condition using a multivariate population decoder.
%   Trials are iteratively withheld and the performance of classifiers trainied 
%   on the remaining N-1 trials is evaluated in successfully determining the 
%   presence/absence of the distractor in the witheld trial.
% - Population information content is then compared between genotypes as a
%   function of the size of randomly drawn pools (ensembles) of neurons
%   used for classification.
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
%% Setup
clear all; close all

if ispc
    pat= 'F:\ModellingResults\PatternSeparation';
else ismac
    pat = '/Volumes/5TB externa/ModellingResults/PatternSeparation/';
end
cd(pat)
stimFreqs=[5 10 20 50];
stimPosition = [0 1 2 3 4];
noModels=5;
noStimConds=length(stimPosition);
noTrials=10;
Ne=800;
% recFolderWT = uigetdir(pat,'Select WT folder');
% recFolderKO = uigetdir(pat,'Select KO folder');
%% Decode WT - Percent correct predictions using Leave-one-out decoding
load([pat filesep 'WTmeta.mat'])
nrep = 10;
% drawSize = Ne/100;
drawSize = [10 50 100 250 500];
for iDraw = 1:length(drawSize)
    for iFreq = 1:length(stimFreqs)
        for iModel=1:noModels
            for iPattern = 2:noStimConds
                FR = [cell2mat(meta.results.SDF{iFreq}{iModel}{1}),cell2mat(meta.results.SDF{iFreq}{iModel}{iPattern})]';
                outcome = [ones(noTrials,1);2*ones(noTrials,1)];
                if drawSize(iDraw)<Ne
                    for iRep=1:nrep
                        fprintf('Decoding Freq %d of %d // model %d of %d // Pattern %d of %d // Repeat no. %d of %d // Drawing %d neurons\n',iFreq,length(stimFreqs),iModel,noModels,iPattern-1,noStimConds-1,iRep,nrep, drawSize(iDraw))
                        rs=randsample(Ne,drawSize(iDraw)); % random draw of units
                        Accuracy(iRep,:) = DecodingError(FR(:,rs),outcome,0.05);
                    end
                else
                    fprintf('Decoding Freq %d of %d // model %d of %d // Pattern %d of %d // Repeat no. 1 of 1 // Drawing all neurons\n',iFreq,length(stimFreqs),iModel,noModels,iPattern-1,noStimConds-1)
                    Accuracy(1,:) = DecodingError(FR(:,1:Ne),outcome,0.05);
                end                    
                Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw}(:,iModel)=1-mean(Accuracy,1); % Collapse results across all draws
            end
        end
    end
end
save([pat filesep 'WTdecodingClassError_UnitNos.mat'],'Decoding','-v7.3')


clear meta Decoding FR
%% Decode KO - Percent correct predictions using Leave-one-out decoding
load([pat filesep 'KOmeta.mat'])
nrep = 10;
% drawSize = Ne/100;
drawSize = [10 50 100 250 500];
for iDraw = 1:length(drawSize)
    for iFreq = 1:length(stimFreqs)
        for iModel=1:noModels
            for iPattern = 2:noStimConds
                FR = [cell2mat(meta.results.SDF{iFreq}{iModel}{1}),cell2mat(meta.results.SDF{iFreq}{iModel}{iPattern})]';
                outcome = [ones(noTrials,1);2*ones(noTrials,1)];
                if drawSize(iDraw)<Ne
                    parfor iRep=1:nrep
                        fprintf('Decoding Freq %d of %d // model %d of %d // Pattern %d of %d // Repeat no. %d of %d // Drawing %d neurons\n',iFreq,length(stimFreqs),iModel,noModels,iPattern-1,noStimConds-1,iRep,nrep, drawSize(iDraw))
                        rs=randsample(Ne,drawSize(iDraw)); % random draw of units
                        Accuracy(iRep,:) = DecodingError(FR(:,rs),outcome,0.05);
                    end
                else
                    fprintf('Decoding Freq %d of %d // model %d of %d // Pattern %d of %d // Repeat no. 1 of 1 // Drawing all neurons\n',iFreq,length(stimFreqs),iModel,noModels,iPattern-1,noStimConds-1)
                    Accuracy(1,:) = DecodingError(FR(:,1:Ne),outcome,0.05);
                end                    
                Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw}(:,iModel)=1-mean(Accuracy,1); % Collapse results across all draws
            end
        end
    end
end
save([pat filesep 'KOdecodingClassError_UnitNos.mat'],'Decoding','-v7.3')
%% calculate average decoding power vs. no units drawn
% Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw}(iModel,:)
drawSize = [10 50 100 250 500];

Decoding.ClassificationError_ = [];
Decoding.ClassificationError_mean = [];
Decoding.ClassificationError_SEM  = [];
for iFreq = 1:length(stimFreqs)
    for iDraw=1:5
        for iPattern = 1:noStimConds-1
             Decoding.ClassificationError_{iFreq}{iDraw}(iPattern,:) = mean(Decoding.ClassificationError{iFreq}{iPattern}{iDraw},2);
        end
        Decoding.ClassificationError_{iFreq}{iDraw}(Decoding.ClassificationError_{iFreq}{iDraw}==0)=0.5;
        Decoding.ClassificationError_mean{iFreq}{iDraw} = nanmean(Decoding.ClassificationError_{iFreq}{iDraw});
        Decoding.ClassificationError_SEM{iFreq}{iDraw}  = nansem(Decoding.ClassificationError_{iFreq}{iDraw});

    end
     Decoding.ClassificationError_mean{iFreq}=cell2mat(Decoding.ClassificationError_mean{iFreq}');
     Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}==0.5)=NaN;
     Decoding.ClassificationError_mean{iFreq}=nanmedian(Decoding.ClassificationError_mean{iFreq},2);
    
     Decoding.ClassificationError_SEM{iFreq}=cell2mat(Decoding.ClassificationError_SEM{iFreq}');
     Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}==0.5)=NaN;
     Decoding.ClassificationError_SEM{iFreq}=nanmean(Decoding.ClassificationError_SEM{iFreq},2);
end
    Decoding.ClassificationError_mean=cell2mat(Decoding.ClassificationError_mean);
    Decoding.ClassificationError_mean(isnan(Decoding.ClassificationError_mean))=0.5;
    Decoding.ClassificationError_SEM=cell2mat(Decoding.ClassificationError_SEM);
    figure; hold on
    for iFreq = 1:length(stimFreqs)
        errorbar(drawSize,Decoding.ClassificationError_mean(:,iFreq),Decoding.ClassificationError_SEM(:,iFreq))
    end

    set(gca,'XScale','log'); 
    grid on; axis square
    axis([min(drawSize) max(drawSize) 0 1])
    xlabel('No. neurons drawn')
    ylabel('Average decoding performance')
    legend(strcat(cellfun(@num2str,num2cell(stimFreqs),'UniformOutput',0),repmat({'Hz stimulation.'},1,numel(stimFreqs))),'Location','northwest'); legend boxoff
%% calculate average decoding power vs. no units drawn - WT and KO
% Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw}(iModel,:)
drawSize = [10 50 100 250 500];
 load([pat filesep 'WTdecodingClassError_UnitNos.mat'])
Decoding.ClassificationError_ = [];
Decoding.ClassificationError_mean = [];
Decoding.ClassificationError_SEM  = [];
for iFreq = 1:length(stimFreqs)
    for iDraw=1:5
        for iPattern = 1:noStimConds-1
            temp = Decoding.ClassificationError{iFreq}{iPattern}{iDraw};
            temp(temp==0.5)=NaN;
            Decoding.ClassificationError_{iFreq}{iDraw}(iPattern,:) = mean(temp,2);
        end
        Decoding.ClassificationError_mean{iFreq}{iDraw} = nanmedian(Decoding.ClassificationError_{iFreq}{iDraw});
        Decoding.ClassificationError_SEM{iFreq}{iDraw}  = nansem(Decoding.ClassificationError_{iFreq}{iDraw});
        
    end
    Decoding.ClassificationError_mean{iFreq}=cell2mat(Decoding.ClassificationError_mean{iFreq}');
    Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}==0.5)=NaN;
    Decoding.ClassificationError_mean{iFreq}=nanmean(Decoding.ClassificationError_mean{iFreq},2);
    
    Decoding.ClassificationError_SEM{iFreq}=cell2mat(Decoding.ClassificationError_SEM{iFreq}');
    Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}==0.5)=NaN;
    Decoding.ClassificationError_SEM{iFreq}=nanmean(Decoding.ClassificationError_SEM{iFreq},2);
end
Decoding.ClassificationError_mean=cell2mat(Decoding.ClassificationError_mean);
Decoding.ClassificationError_mean(isnan(Decoding.ClassificationError_mean))=0.5;
Decoding.ClassificationError_SEM=cell2mat(Decoding.ClassificationError_SEM);
WTDecoding =Decoding;

load([pat filesep 'KOdecodingClassError_UnitNos.mat'])
Decoding.ClassificationError_ = [];
Decoding.ClassificationError_mean = [];
Decoding.ClassificationError_SEM  = [];
for iFreq = 1:length(stimFreqs)
    for iDraw=1:5
        for iPattern = 1:noStimConds-1
            temp = Decoding.ClassificationError{iFreq}{iPattern}{iDraw};
            temp(temp==0)=NaN;
            Decoding.ClassificationError_{iFreq}{iDraw}(iPattern,:) = nanmean(temp,2);
        end

        Decoding.ClassificationError_mean{iFreq}{iDraw} = nanmedian(Decoding.ClassificationError_{iFreq}{iDraw});
        Decoding.ClassificationError_SEM{iFreq}{iDraw}  = nansem(Decoding.ClassificationError_{iFreq}{iDraw});
        
    end
    Decoding.ClassificationError_mean{iFreq}=cell2mat(Decoding.ClassificationError_mean{iFreq}');
    Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}==0.5)=NaN;
    Decoding.ClassificationError_mean{iFreq}=nanmean(Decoding.ClassificationError_mean{iFreq},2);
    
    Decoding.ClassificationError_SEM{iFreq}=cell2mat(Decoding.ClassificationError_SEM{iFreq}');
    Decoding.ClassificationError_mean{iFreq}(Decoding.ClassificationError_mean{iFreq}<0.5)=NaN;
    Decoding.ClassificationError_SEM{iFreq}=nanmean(Decoding.ClassificationError_SEM{iFreq},2);
end
Decoding.ClassificationError_mean=cell2mat(Decoding.ClassificationError_mean);
Decoding.ClassificationError_mean(isnan(Decoding.ClassificationError_mean))=0.5;
Decoding.ClassificationError_SEM=cell2mat(Decoding.ClassificationError_SEM);
KODecoding =Decoding;

    figure; 
    
    for iFreq = 3:length(stimFreqs)
        subplot(1,length(stimFreqs),iFreq);hold on
%         subaxis(1,2,iFreq-2);hold on
        errorbar(drawSize,100*WTDecoding.ClassificationError_mean(:,iFreq),100*Decoding.ClassificationError_SEM(:,iFreq),'b','LineWidth',1.5)
        errorbar(drawSize,100*KODecoding.ClassificationError_mean(:,iFreq),100*Decoding.ClassificationError_SEM(:,iFreq),'r','LineWidth',1.5)
        plot([min(drawSize) max(drawSize)], [50 50],':k','LineWidth',1.5);
        set(gca,'XScale','log'); 
        grid on; axis square; box on
%     axis([min(drawSize) max(drawSize)+2 50 100])
    axis([0 550 35 100])
    set(gca,'Xtick',drawSize,'XtickLabel',drawSize,'XTickLabelRotation',45)
    xlabel('Population size (no. neurons)')
    if iFreq==3
            ylabel({'Mean population performance';'(% Decoders Correct)'})
    else
        set(gca,'YTicklabel',[])

    end
    end
%% Plot group decoding for each distractor (time resolved)
%     legend(strcat(cellfun(@num2str,num2cell(stimFreqs),'UniformOutput',0),repmat({'Hz stim'},1,numel(stimFreqs))),'Location','northwest'); legend boxoff
x_temp = 1:size(WT.Decoding.ClassificationError{1}{1},1);

StimTimes{1} = [10 210 410 610 810];
StimTimes{2} = [10 110 210 310 410];
StimTimes{3} = [10 60  110 160 210];
StimTimes{4} = [10 30  50  70  90];

Distractors{1} = [110 210 510 710];
Distractors{2} = [60 160 260 360];
Distractors{3} = [35 85 135 185];
Distractors{4} = [20 40 50 80];
xlims_ = [1000 500 400 200];
offset = 1;

WT = load([pat filesep 'WTdecodingClassError.mat']);
KO = load([pat filesep 'KOdecodingClassError.mat']);

figure('color','w');
for iFreq = 1:length(stimFreqs)
    subplot(1,length(stimFreqs),iFreq); hold on
    for iPattern = 1:noStimConds-1

        
        
        
        
        if iFreq>2
            plot(offset*iPattern+mean(WT.Decoding.ClassificationError{iFreq}{iPattern},2),'b','LineWidth',1.5)
            ciplot(offset*iPattern+mean(WT.Decoding.ClassificationError{iFreq}{iPattern},2)+nansem(WT.Decoding.ClassificationError{iFreq}{iPattern},2),...
                    offset*iPattern+mean(WT.Decoding.ClassificationError{iFreq}{iPattern},2)-nansem(WT.Decoding.ClassificationError{iFreq}{iPattern},2),...
            x_temp,'b')
        end
            plot(offset*iPattern+mean(KO.Decoding.ClassificationError{iFreq}{iPattern},2),'r','LineWidth',1.5)
             ciplot(offset*iPattern+mean(KO.Decoding.ClassificationError{iFreq}{iPattern},2)+nansem(KO.Decoding.ClassificationError{iFreq}{iPattern},2),...
                    offset*iPattern+mean(KO.Decoding.ClassificationError{iFreq}{iPattern},2)-nansem(KO.Decoding.ClassificationError{iFreq}{iPattern},2),...
            x_temp,'r')
%         plot(WT.Decoding.Ft2ciH_mean{iFreq}(iPattern,:),'b')
%         plot(KO.Decoding.Ft2ciH_mean{iFreq}(iPattern,:),'r')

        plot([Distractors{iFreq}(iPattern); Distractors{iFreq}(iPattern)],...
        [offset*iPattern*ones(1,length(StimTimes{iFreq}));offset*iPattern*ones(1,length(StimTimes{iFreq}))],'color',[0.9 0.6 0.3],'LineWidth',4)%[0.9 0.6 0.3]
        plot([StimTimes{iFreq}; StimTimes{iFreq}],[offset*iPattern-offset;offset*iPattern-offset],'color','k','LineWidth',4)%[0.9 0.6 0.3]
    end
    if iFreq>2
        axis([0 xlims_(iFreq)  0 offset*5])
    else
        axis([0 xlims_(iFreq) 0 offset*5])
    end
%     axis off
end
% 
% ylabel('Prediction accuracy (F-score)')
% xlabel('Time (ms)')
%% Population analysis
WTDecoding = load('WTdecodingClassError_UnitNos.mat');
KODecoding = load('KOdecodingClassError_UnitNos.mat');
%%

x_temp = 1:size(WTDecoding.Decoding.ClassificationError{1}{1}{1},1);

StimTimes{1} = [10 210 410 610 810];
StimTimes{2} = [10 110 210 310 410];
StimTimes{3} = [10 60  110 160 210];
StimTimes{4} = [10 30  50  70  90];

Distractors{1} = [110 210 510 710];
Distractors{2} = [60 160 260 360];
Distractors{3} = [35 85 135 185];
Distractors{4} = [20 40 50 80];

iFreq = 3;
iDraw = 1;

NBS = 200;


figure('color','w')
for iPattern=2:5
    subplot(4,1,iPattern-1) ; hold on
    WT_ = WTDecoding.Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw};
    KO_=  KODecoding.Decoding.ClassificationError{iFreq}{iPattern-1}{iDraw};
    
    WT_(1:Distractors{iFreq}(iPattern-1),:)=0.5;
    KO_(1:Distractors{iFreq}(iPattern-1),:)=0.5;
    
    [p,~] = permtest2vec(WT_,KO_,NBS);
        %     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
        a = nan(size(p));a(p)=1;
        a(mean(WT_,2)==mean(KO_,2))=NaN;
        plot(1:1000,1.3*a,'-k','LineWidth',2.5)
ciplot(nanmean(WT_,2)+nanstd(WT_,2),nanmean(WT_,2)-nanstd(WT_,2),1:1000,'b',0.7)
ciplot(nanmean(KO_,2)+nanstd(KO_,2),nanmean(KO_,2)-nanstd(KO_,2),1:1000,'r',0.7)   
plot(1:1000,mean(WT_,2),'b','LineWidth',1.5)
plot(1:1000,mean(KO_,2),'r','LineWidth',1.5)

plot([StimTimes{iFreq}; StimTimes{iFreq}],[0.3*ones(length({StimTimes}));0.4*ones(length({StimTimes}))],'k','LineWidth',4)
plot([Distractors{iFreq}(iPattern-1); Distractors{iFreq}(iPattern-1)],[0.3;0.4],'color',[0.9 0.6 0.3],'LineWidth',4)
% axis([0 200 0.25 1.5]);% 50Hz
axis([0 450 0.25 1.5]);% 200Hz
axis off
end



