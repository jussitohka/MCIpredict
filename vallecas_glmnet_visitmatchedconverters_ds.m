% The main glmnet script to do the data-analysis
%    Copyright (C) Jussi Tohka, 2020
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along
%    with this program; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

clear
close all
if isunix
   datadir = '/research/work/justoh/vallecas_study/data';
    matlabdir = '/research/users/justoh/matlab'; 
    savedir = '/research/work/justoh/vallecas_study/results';
else
    datadir = 'C:\Users\justoh\Data\vallecas_study';
    matlabdir = 'C:\Users\justoh\matlab';
    savedir = 'C:\Users\justoh\Results\vallecas_study';
end
addpath(fullfile(matlabdir,'glmnet_matlab2'));
% addpath(fullfile(matlabdir,'libsvm-3.21','matlab'));
addpath(fullfile(matlabdir,'export_fig-master'));

ctrlidxfn =  'matched_controlidx_newconverters';
load(fullfile(datadir,'vdata_visitmatched_v4'))
load(fullfile(datadir,'vdata_visitmatched_v4_varnames'));
lambdavec = 0.01:0.01:0.5;
alphavec = 0:0.01:1;


family = 'binomial';
type = 'auc';
opt = glmnetSet;
opt.lambda = lambdavec;

feature_type = {'Apoe+Demographics','Demographics','Neuropsychology','MRI','Neuropsych+demographics','All'};
nfeature_types = length(feature_type);


% create the complete data set
vdata.data{5} = [vdata.data{2} vdata.data{3}];
vdata.data{6} = [vdata.data{1} vdata.data{3} vdata.data{4}];

varnames{5} = [varnames{2} varnames{3}];
varnames{6} = [varnames{1} varnames{3} varnames{4}];

% standardize the whole data 
idx12 = find(vdata.label < 3);
% standardize the whole data 
idx12 = find(vdata.label < 3);
for k = 1:nfeature_types
    [~,mu{k},sigma{k}] = zscore(vdata.data{k}(idx12,:));
    tmp = bsxfun(@minus,vdata.data{k},mu{k});
    vdata.data{k} = bsxfun(@rdivide,tmp,sigma{k});
end



y = vdata.label(idx12);
weights = zeros(size(y));
weights(y == 1) = 1/sum(y == 1);% weighting the data
weights(y == 2) = 1/sum(y == 2);
weights = weights*length(y)/2;  
opt.weights = weights;

% either run the glmnet or if it has been already run load the results
if isunix | ~exist(fullfile(savedir,'results_visitmatchedconverters_withpreds.mat'),'file') 
    for k = 1:nfeature_types
        x = vdata.data{k}(idx12,:);
        for a = 1:length(alphavec)
            disp([k a])
            opt.alpha = alphavec(a);
            cverra{k,a} = cvglmnet_jt(x,y,family,opt,type,10,vdata.foldid12,0,1,0,25);
         %   fs(:,i) = glmnetPredict(cverra.glmnet_fit,x,cverra.lambda_min,'nonzero');
        end
    end
    save(fullfile(savedir,'results_visitmatchedconverters_withpreds'),'cverra');
else 
    load(fullfile(savedir,'results_visitmatchedconverters_withpreds.mat'))
end

y = y == 2;

for k = 1:length(varnames)
    varnames{k} = [ {'Const'}, varnames{k}];
end

for k = 1:nfeature_types
    for a = 1:length(alphavec)
         [bestmodel{k}(a), lambdaidx{k}(a)] = max(cverra{k,a}.cvm);
         bestlambda{k}(a) = cverra{k,a}.glmnet_fit.lambda(lambdaidx{k}(a)); 
    end
end

for k = 1:nfeature_types
    [maxauc(k),aucidx(k)] = max(bestmodel{k});
    % find aucs of the best model during repeated CV
    aucs_repeats = cverra{k,aucidx(k)}.cvraw(:,lambdaidx{k}(aucidx(k)));
    % find the median value and its index
    % function median does not give the index
    [cvs,idx]  = sort(aucs_repeats);
    medianidx(k) = idx(ceil(length(aucs_repeats)/2));
    median_auc(k) = cvs(ceil(length(aucs_repeats)/2));
    valpred{k} = cverra{k,aucidx(k)}.fit_preval{medianidx(k)}(:,lambdaidx{k}(aucidx(k)));
end

% cverra = rmfield(cverra,'fit_preval');
% save(fullfile(savedir,'results_visitmatchedconverters_withpreds_small'),'cverra','valpred');

idx34 = find(vdata.label > 2 & vdata.label < 5);
ytest = (vdata.label(idx34) - 2) > 1.5;
% classes are 1 and 2, but note that pred{k} is the class probability so below 
% conversion gives classes 0 and 1  
% class 0 -> converters (positive)
% class 1 -> controls (negative)
% this is equal to running glmnetPredict with 'class'
poslabel = 0;
for k = 1:nfeature_types
   xtest = vdata.data{k}(idx34,:);
   
   bestobject{k} = glmnetPredict(cverra{k,aucidx(k)}.glmnet_fit,[],bestlambda{k}(aucidx(k)),'coefficients', [],[]);
   pred{k} = glmnetPredict(cverra{k,aucidx(k)}.glmnet_fit,xtest,bestlambda{k}(aucidx(k)),'response', [],[]);
   predlabels{k} = pred{k} > 0.5;
   coefidx{k} = find(bestobject{k});
   coefvalues{k} = bestobject{k}(coefidx{k});
end

makeplot = 0;
boot = 0;
cv = 50; % this should be 0 or 50, 50 for split-half calibration
nfolds = 2;
foldfn = strcat('vallecas_test_folds',num2str(nfolds));
rc = {'r','g','b','b:','m','y','k'};
if cv
    if ~exist(fullfile(savedir,foldfn),'file')
        for iter = 1:cv
            foldid{iter} = balanced_crossval(ytest,nfolds,[],0,0);
        end
        save(fullfile(savedir,foldfn),'foldid');
    else
        load(fullfile(savedir,foldfn))
    end
end   
for k = 1:nfeature_types
   cm{k} = confusionmat(ytest,predlabels{k});
   [sen(k),spec(k),acc(k)] = senspec(ytest,predlabels{k},poslabel);
   valpredlabels{k} = valpred{k} > 0.5;
   [valsen(k),valspec(k),valacc(k)] = senspec(y,valpredlabels{k},poslabel);
   [rocX{k},rocY{k},T{k},testauc{k},bestopt] = perfcurve(ytest,1 - pred{k},0); %,'nboot',5000,'xvals','all');
   if cv
       for iter = 1:cv
           
           for f = 1:nfolds
               which = foldid{iter} == f;
               [cvrocX,cvrocY,cvT,cvtestauc,cvbestopt] = perfcurve(ytest(~which),1 - pred{k}(~which),0);
               cvbest_thr(k,iter,f) = cvT((cvrocX==cvbestopt(1))&(cvrocY==cvbestopt(2)));
               cvoptpredlabels = pred{k}(which) > (1 - cvbest_thr(k,iter,f));
               [cvoptsen(k,iter,f),cvoptspec(k,iter,f),cvoptacc(k,iter,f)] = senspec(ytest(which),cvoptpredlabels,poslabel);
           end
       end
       tmp = squeeze(cvbest_thr(k,:,:));
       median_thr(k) = median(tmp(:));
   end
   
   [valrocX{k},valrocY{k},Tv,valauc{k}] = perfcurve(y,1 - valpred{k},0,'Weights',weights); % ,'nboot',5000,'xvals','all');
   best_thr(k) = T{k}((rocX{k}==bestopt(1))&(rocY{k}==bestopt(2)));
   optpredlabels{k} = pred{k} > (1 - best_thr(k));
   [optsen(k),optspec(k),optacc(k)] = senspec(ytest,optpredlabels{k},poslabel);
   if boot
       [brocX{k},brocY{k},T{k},testauc{k}] = perfcurve(ytest,1 - pred{k},0,'nboot',5000,'xvals','all');
       
   end
  
end


if makeplot   % the fonts and colors have been changed in the paper manually  
    figure
    hold
    for k = 1:nfeature_types
       
       plot(rocX{k},rocY{k},rc{k},'LineWidth',2)
       
    end
    plot([0 1],[0 1],'c','LineWidth',2)
    legend([feature_type([1:3 5 7]),'Constant'],'Location','SouthEast','FontSize',12)
    xlabel('FPR','FontSize',12)
    ylabel('TPR','FontSize',12)
    export_fig(fullfile(savedir,strcat('testroc_1year_42conv_42controls_v2','.png')))
    saveas(gcf,fullfile(savedir,'testroc_1year_42conv_42controls_v2'),'fig')
    close
    
    for k = 1:nfeature_types
            figure
            hold
            plot(rocX{k},rocY{k},rc{k},'LineWidth',2)
            plot(valrocX{k},valrocY{k},[rc{k},':'],'LineWidth',2)
            plot([0 1],[0 1],'c','LineWidth',2)
            title(feature_type{k})
            legend({'Test ROC','Validation ROC','Constant'},'Location','SouthEast','FontSize',12)
            export_fig(fullfile(savedir,strcat('valtestroc_1year_',feature_type{k},'.png')))
            saveas(gcf,fullfile(savedir,strcat('valtestroc_1year',feature_type{k})),'fig');
      
    end
end

for k = 1:nfeature_types
    disp('************************************************');
    disp(feature_type{k});
    disp('************************************************');
    disp('Model: Std value, real value ');
    % zi = (xi - mi)/si;
    % y = sum(ai*zi) + c;
    % y = =sum(ai*xi/si) - sum(ai*mi/si) + c;
    for i = 1:length(coefidx{k})
        if i > 1
            rv = coefvalues{k}(i)/sigma{k}(coefidx{k}(i) - 1); % -1 is because the constant term is 
                                                                  % added to the beginning of coefficients
                                                                  % but not to sigma and mu                                   
            disp([varnames{k}{coefidx{k}(i)}, ' : ', num2str(coefvalues{k}(i)),' ',num2str(rv)]);
        else
            rv = (-1)*sum(mu{k}(coefidx{k}(2:end) - 1).*(coefvalues{k}(2:end)')./sigma{k}(coefidx{k}(2:end) - 1)) + coefvalues{k}(1);
            disp([varnames{k}{coefidx{k}(i)}, ' : ', num2str(coefvalues{k}(i)),' ',num2str(rv)]);
        end
    end 
    disp(['Best Cross-validated AUC:', num2str(maxauc(k)), ';Alpha:',num2str(alphavec(aucidx(k))),';Lambda:',num2str(bestlambda{k}(aucidx(k)))]);
    disp(['Test AUC + 95% CI:',num2str(testauc{k})]);
    disp(['Test sensitivity:',num2str(sen(k))])
    disp(['Test accuracy:',num2str(acc(k)),' Specificity:',num2str(spec(k))])
    disp(['Optimal acc', num2str(optacc(k)), ' Sen:' num2str(optsen(k)), ' Spec:', num2str(optspec(k))]);
    disp(['Optimal threshold:' num2str(best_thr(k))]);
    disp('Sensitivity is % of correctly identified converters.')
    disp('Specificity is % of correctly identified controls.')
end


% Prepare data for star http://melolab.org/star/roc_analysis.php

posidx = find(ytest == 0);
negidx = find(ytest == 1);
pospred = [];
negpred = [];
for k = 1:nfeature_types
    pospred =[pospred, 1 - pred{k}(posidx)];
    negpred = [negpred, 1 - pred{k}(negidx)];
end






