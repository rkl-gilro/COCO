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
Ntrials=25:50:125;
Nstimuli= [10 15 20];
ntest=2;
testValues=linspace(.2,.8,ntest);

% testValue=.62;
sigma=.2;
Sigmas=linspace(.1,.5,2);
%%% SIMULATE
nreps=1;
c=0;
redo=1;
if redo ==1
for rep=1:nreps
    for sg=1:length(Sigmas)
        sigma=Sigmas(sg);
        for ns=1:length(Nstimuli)
            nstimuli=Nstimuli(ns);
            CompetitorsValues= Scale(sort(rand(nstimuli,1)));
      
           
            competitorIndices=nchoosek(1:nstimuli,5);
            competitorIndicesPairs=nchoosek(1:nstimuli,2);
            


            for ntr=1:length(Ntrials)
                ntrials=Ntrials(ntr);
                
                error pd
                
                
    for test=1:ntest
                    testValue=testValues(test);
                    
                    
                    responses=zeros(size(competitorIndices,1),1);
                    numbertrials=zeros(size(competitorIndices,1),1);
                    for trial=1:ntrials
                       
                        compgroupind= ceil(rand*size(competitorIndices,1) + 0.000001);
                        inds=competitorIndices(compgroupind,:);
                        s= CompetitorsValues(competitorIndices(compgroupind,:));
                        ds= abs(s'-testValue+randn(1,length(s))*sigma);
%                         d1= abs(s1-testValue+randn*sigma);
%                         d2= abs(s2-testValue+randn*sigma);

winnerInd= find(ds==min(ds));
LoosersInd= find(ds~=min(ds));
pairs=[repmat(winnerInd,length(LoosersInd),1) LoosersInd'];
%                         resp= d1<d2;
for pp=1:size(pairs,1)
   pos=find(sum( repmat(pairs(pp,:),size(competitorIndicesPairs,1),1)-competitorIndicesPairs,2)==0);
% include also he opposite pairing order, so that looses are also counted
   responses(pos)=responses(pos)+resp;
                        numbertrials(pair)=numbertrials(pair)+1;
   
end

                        
                    end
% error pd
                    fl=['T' num2str(test) 'TR' num2str(ntr) 'NS' num2str(ns) 'SG' num2str(sg) 'NR', num2str(rep)];
                 saveMLDSColorSelectionMT(competitorIndices,responses,numbertrials, nstimuli,fl);
   
%                     [targetCompetitorFit, logLikelyFit, predictedResponses] = MLDSColorSelection(competitorIndices,responses,numbertrials, nstimuli);
                    %   figure
%                     MLDSColorSelectionPlot(competitorIndices,responses,numbertrials,targetCompetitorFit,predictedResponses, 0);
   end
            for test=1:ntest
 fl=['T' num2str(test) 'TR' num2str(ntr) 'NS' num2str(ns) 'SG' num2str(sg) 'NR', num2str(rep)];

load(fl)
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
% save SIMULATIONS_RANDOM correlationEstimates
else % read saved files



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
% ylim([0 1])
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
tmp(:,test)=[correlateEstimates{test,:,ns,sg,rep}];
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
% ylim([.8 1])
 end
