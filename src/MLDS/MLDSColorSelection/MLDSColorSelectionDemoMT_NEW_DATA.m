% MLDSColorSelectionDemo.m
%
% Demonstrates MLDS Fitting procedure for a data set.
% Initially used as a test bed for improving search algorithm.
% Currently used as demo and also to debug the algorithm for
% any "problematic data set".
%
% The work is done by other routines in this folder.
%
% Requires optimization toolbox.
%
% 5/3/12  dhb  Wrote it. Added scatterplot of predicted versus measured probs added.
% 6/13/13  ar  Clean up, added comments

% Initialize
clear ; close all;


% % % %     simulate responses
rng(1)
Ntrials=50:50:150;
Nstimuli= 10:4:22;
ntest=5;
testValues=linspace(.2,.6,ntest);

% testValue=.62;
sigma=.2;
Sigmas=linspace(.1,.5,3);
%%% SIMULATE
nreps=1;
c=0;
redo=0;
if redo ==1
for rep=1:nreps
    for sg=1:length(Sigmas)
        sigma=Sigmas(sg);
        for ns=1:length(Nstimuli)
            nstimuli=Nstimuli(ns)
            CompetitorsValues= Scale(sort(rand(nstimuli,1)));
            competitorIndices=nchoosek(1:nstimuli,2);
            
            for ntr=1:length(Ntrials)
                ntrials=Ntrials(ntr);
                
                for test=1:ntest
                    testValue=testValues(test);
                    responses=zeros(size(competitorIndices,1),1);
                    numbertrials=zeros(size(competitorIndices,1),1);
                    for trial=1:ntrials
                        pair= ceil(rand*size(competitorIndices,1) + 0.000001);
                        s1= CompetitorsValues(competitorIndices(pair,1));
                        s2= CompetitorsValues(competitorIndices(pair,2));
                        d1= abs(s1-testValue+randn*sigma);
                        d2= abs(s2-testValue+randn*sigma);
                        resp= d1<d2;
                        responses(pair)=responses(pair)+resp;
                        numbertrials(pair)=numbertrials(pair)+1;
                    end
                    
                    
                    [targetCompetitorFit, logLikelyFit, predictedResponses] = MLDSColorSelection(competitorIndices,responses,numbertrials, nstimuli);
                    %   figure
%                     MLDSColorSelectionPlot(competitorIndices,responses,numbertrials,targetCompetitorFit,predictedResponses, 0);
                    targetCompetitorFit_z=zscore(targetCompetitorFit);
                    estimates(test)=targetCompetitorFit_z(1);
                    correlateEstimates{test,ntr,ns,sg,rep}= corr(targetCompetitorFit(2:end)', CompetitorsValues);
                    correlateResponses{test,ntr,ns,sg,rep}= corr(predictedResponses', responses./numbertrials);
                c=c+1;
pc = c/(ntest*length(Ntrials)*length(Nstimuli)*length(Sigmas)*nreps);
pc = round(pc*100,2);
disp([num2str(pc) '%'])

end
                correlationEstimates{ntr,ns,sg,rep}=corr(estimates(:),testValues(:));
            end
            
        end
    end
end
% figure
% plot(estimates,testValues,'ko')
save SIMULATIONS_RANDOM correlationEstimates
else
load SIMULATIONS_RANDOM
end

symbs={'-','--',':'};
 for sg=1:length(Sigmas)
subplot(1,length(Sigmas),sg)
    for rep=1:nreps
   
        for ns=1:length(Nstimuli)
%         for ntr=1:length(Ntrials)
%             col=ones(1,3)*(ntr/length(Ntrials))*.8;
col=ones(1,3)*(ns/length(Nstimuli))*.8;
%             plot( Nstimuli, [  correlationEstimates{ntr,:,sg,rep}],['ko' symbs{sg}],'markerfacecolor',col,'markeredgecolor',col)
%   plot( Ntrials, [  correlationEstimates{:,ns,sg,rep}],['ko' symbs{sg}],'markerfacecolor',col,'markeredgecolor',col,'linewidth',2,'markersize',14)

plot( Ntrials, [  correlationEstimates{:,ns,sg,rep}],'ko-','markerfacecolor',col,'markeredgecolor',col,'linewidth',2,'markersize',14)
                  
 hold on
        end
    end
   
% xlabel('#Competitors','FontSize',15)
xlabel('#Trials','FontSize',15)

ylabel('Correlation (estimated ~ real target)','FontSize',15)
set(gca,'FontSize',15,'LineWidth',3)
box off 
axis square
xlim([40 170])
ylim([0 1])
 end


figure
 for sg=1:length(Sigmas)
subplot(1,length(Sigmas),sg)
    for rep=1:nreps
   
        for ns=1:length(Nstimuli)
%         for ntr=1:length(Ntrials)
%             col=ones(1,3)*(ntr/length(Ntrials))*.8;
col=ones(1,3)*(ns/length(Nstimuli))*.8;
%             plot( Nstimuli, [  correlationEstimates{ntr,:,sg,rep}],['ko' symbs{sg}],'markerfacecolor',col,'markeredgecolor',col)
%   plot( Ntrials, [  correlationEstimates{:,ns,sg,rep}],['ko' symbs{sg}],'markerfacecolor',col,'markeredgecolor',col,'linewidth',2,'markersize',14)
  
for test=1:ntest  
tmp(:,test)=[correlationEstimates{test,:,ns,sg,rep}];
end
tmp=mean(tmp,2);
  
plot( Ntrials, tmp,'ko-','markerfacecolor',col,'markeredgecolor',col,'linewidth',2,'markersize',14)
                  
 hold on
        end
    end
   
% xlabel('#Competitors','FontSize',15)
xlabel('#Trials','FontSize',15)

ylabel('Correlation (estimated ~ real competitors)','FontSize',15)
set(gca,'FontSize',15,'LineWidth',3)
box off 
axis square
xlim([40 170])
ylim([.8 1])
 end
