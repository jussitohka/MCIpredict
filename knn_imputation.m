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

% Missing data imputation by k nearest neighbours
% data is nsubj times nvar data matrix
% subjidx are the indexes of subjects that can be used for imputation (aka
% training set)
% varidx are the indexes of variables that cannot be used for imputation
% k is the number of nearest neighbours
% Reference:
% Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, 
% David Botstein, Russ B. Altman, Missing value estimation methods for DNA microarrays , 
%Bioinformatics, Volume 17, Issue 6, June 2001, Pages 520–525

function data = knn_imputation(data,subjidx,varidx,k)
% first find the subjects with missing data
[iii,jjj] = find(isnan(data));
uiii = unique(iii); % unique list of subjects with missing values
% remove these subjs from subjidx
subjidx = setdiff(subjidx,iii);

% compute standradized data for the included subjects
[zdata,mu,sigma] = zscore(data(subjidx,:));
zdata_complete = bsxfun(@rdivide,(data - mu),sigma);
data_rest = data(subjidx,:);
for i = 1:length(uiii)
    ujjj = find(isnan(data(uiii(i),:)));
    varidx2 = setdiff(varidx,ujjj);
    idx = knnsearch(zdata_complete(subjidx,varidx2),zdata_complete(uiii(i),varidx2),'K',k);
    for j = 1:length(ujjj)
        data(uiii(i),ujjj(j)) = mean(data_rest(:,ujjj(j)));
    end
end

