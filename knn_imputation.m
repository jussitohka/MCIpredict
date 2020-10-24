% 
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

% Missing data imputation by k nearest neighbours or mean
% requires statistics and machine learning toolbox (nanmean,nanstd)
% data is nsubj times nvar data matrix
% subjidx are the indexes of subjects that can be used for imputation (aka
% training set)
% varidx are the indexes of variables that can be used for imputation
% k is the number of nearest neighbours, give 0 if mean
% Reference:
% Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, 
% David Botstein, Russ B. Altman, Missing value estimation methods for DNA microarrays , 
% Bioinformatics, Volume 17, Issue 6, June 2001, Pages 520–525

function data = knn_imputation(data,subjidx,varidx,k)

subjidx_orig = subjidx; 
% first find the subjects and variables with missing data
[iii,jjj] = find(isnan(data));
% compute how many subjcts are missing data from each variable 
for i = 1:length(varidx)
    nanvars(i)= length(find(isnan(data(:,i))));
end
badvars = find(nanvars > length(subjidx)*0.1);
% remove variables with too many missing subjects from the dataset given to
% knn
varidxb = setdiff(varidx,badvars);

uiii = unique(iii); % unique list of subjects with missing values
% remove these subjs from subjidx
subjidx = setdiff(subjidx,iii);

% construct subjidx for knn
[biii,bjjj] = find(isnan(data(:,varidxb)));
bsubjidx = setdiff(subjidx_orig,unique(biii));
% compute standradized data for the included subjects
% [zdata,mu,sigma] = zscore(data(subjidx,:));
mu = nanmean(data);
sigma = nanstd(data);
zdata_complete = bsxfun(@rdivide,(data - mu),sigma);
data_rest = data(bsubjidx,:);
for i = 1:length(uiii)
    ujjj = find(isnan(data(uiii(i),:)));
    varidx2 = setdiff(varidxb,ujjj);
    bsubjidx2 = setdiff(bsubjidx,uiii(i)); 
    if k > 0
        idx = knnsearch(zdata_complete(bsubjidx2,varidx2),zdata_complete(uiii(i),varidx2),'K',k);
        for j = 1:length(ujjj)
            % nanmean is needed for those variables which are 'badvars'
            
            data(uiii(i),ujjj(j)) = nanmean(data(bsubjidx2(idx),ujjj(j)));
            dbj = data(uiii(i),ujjj(j));
            % if everything else fails take nanmean 
            if isnan(dbj)
                data(uiii(i),ujjj(j)) = nanmean(data(:,ujjj(j)));
            end
        end
    else
        for j = 1:length(ujjj)
            dbj = data(uiii(i),ujjj(j));
            data(uiii(i),ujjj(j)) = nanmean(data(:,ujjj(j)));
        end
    end
end

