% CHANGED BY JUSSI TOHKA 22.3.2020 BY ALLOWING AUC COMPUTATION BY POOLING
function result = cvlognet(object, lambda, x, y, weights, offset, foldid, ...
    type, grouped, keep,repeats)

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('mse','Mean-Squared Error','mae','Mean Absolute Error',...
    'deviance','Binomial Deviance','auc','AUC','class','Misclassification Error');

if strcmp(type,'default')
    type = 'deviance';
end

if ~any(strcmp(type,{'mse','mae','deviance','auc','class'}))
    warning('Only ''deviance'', ''class'', ''auc'', ''mse'' or ''mae''  available for binomial models; ''deviance'' used');
    type = 'deviance';
end

prob_min = 1e-5; prob_max = 1 - prob_min;
nc = size(y);
if nc(2) == 1
    [classes,~,sy] = unique(y);
    nc = length(classes);
    indexes = eye(nc);
    y = indexes(sy,:);
end
N = size(y,1);
nfolds = max(foldid{1});
if (N/nfolds < 10) && strcmp(type,'auc')
    warning(strcat('Too few (< 10) observations per fold for type.measure=''auc'' in cv.lognet; ',...
        'changed to type.measure=''deviance''. Alternatively, use smaller value for nfolds'));
    type = 'deviance';
end
if (N/nfolds < 3) && grouped
    warning(strcat('Option grouped=FALSE enforced in cv.glmnet, ',...
        'since < 3 observations per fold'));
    grouped = false;
end
is_offset = ~isempty(offset);
for k = 1:repeats
     predmat{k} = NaN(size(y,1),length(lambda));
end
nlams = zeros(repeats,nfolds);

for k = 1:repeats
    for i = 1:nfolds
        which = foldid{k}==i;
        fitobj = object{k,i};
        if (is_offset)
            off_sub = offset(which,:);
        else
            off_sub = [];  %a bit different from that in R
        end
        preds{k} = glmnetPredict(fitobj,x(which,:),[],'response',[],off_sub);
        nlami = length(object{k,i}.lambda);
        predmat{k}(which,1:nlami) = preds{k};
        nlams(k,i) = nlami;
    end
end

if strcmp(type,'auc')
    cvraw = NaN(nfolds, length(lambda));
    good = zeros(nfolds, length(lambda));
    for k = 1:repeats
        for i = 1:nfolds
            good(i,1:nlams(k,i)) = 1;
            which = foldid{k} == i;
        end
        for j = 1:nlams(i)
     %           cvraw(i,j) = auc_mat(y(which,:), predmat(which,j), weights(which));
            cvraw(k,j) = auc_mat(y, predmat{k}(:,j), weights);
        end
       
    end
    N = sum(good(1,:),1);
    sweights = zeros(repeats, 1);
    for k = 1:repeats
        sweights(k) = 1;
       % for i = 1:nfolds
       %     sweights(i) = sum(weights(foldid{k}==i));
       % end
    end
    weights = sweights;
else
    ywt = sum(y, 2);
    y = y ./ repmat(ywt,1,size(y,2));
    weights = weights .* ywt;
    N = size(y,1) - sum(isnan(predmat),1);
    yy1 = repmat(y(:,1),1,length(lambda));
    yy2 = repmat(y(:,2),1,length(lambda));
    switch type
        case 'mse'
            cvraw = (yy1 - (1 - predmat)).^2 + (yy2 - (1 - predmat)).^2;
        case 'mae'
            cvraw = abs(yy1 - (1 - predmat)) + abs(yy2 - (1 - predmat));
        case 'deviance'
            predmat = min(max(predmat,prob_min),prob_max);
            lp = yy1.*log(1-predmat) + yy2.*log(predmat);
            ly = log(y);
            ly(y == 0) = 0;
            ly = (y.*ly) * [1;1];
            cvraw = 2 * (repmat(ly,1,length(lambda)) - lp);
        case 'class'
            cvraw = yy1.*(predmat > 0.5) + yy2.*(predmat <= 0.5);
    end
    if (grouped)
        cvob = cvcompute(cvraw, weights, foldid, nlams);
        cvraw = cvob.cvraw;
        weights = cvob.weights;
        N = cvob.N;
    end
end
cvm = wtmean(cvraw,weights);
sqccv = (bsxfun(@minus,cvraw,cvm)).^2;
cvsd = sqrt(wtmean(sqccv,weights)./(N-1));
result.cvm = cvm; result.cvsd = cvsd; result.name = typenames.(type);
result.cvraw = cvraw;
if (keep)
    result.fit_preval = predmat;
end


function result = auc_mat(y, prob, weights)

if nargin < 3 || isempty(weights)
    weights = ones(size(y,1),1);
end
Weights = bsxfun(@times,weights,y);
Weights = Weights(:)';
ny = size(y,1);
Y = [zeros(ny,1);ones(ny,1)];
Prob = [prob; prob];
result = auc(Y, Prob, Weights);


function result = auc(y, prob, w)
mindiff = max(prob) - min(prob);
if (mindiff < eps) | (length(unique(prob)) < 11) % this is tuned for 10-fold CV
    result = 0.5;
else
if isempty(w)
    mindiff = min(diff(unique(prob)));
    pert = unifrnd(0,mindiff/3,size(prob));
    [~,~,rprob] = unique(prob+pert);
    n1 = sum(y); n0 = length(y) - n1;
    u = sum(rprob(y == 1)) - n1*(n1+1)/2;
    result = u / (n1*n0);
else
    [~,op] = sort(prob);
    y = y(op); w = w(op);
    cw = cumsum(w);
    w1 = w(y == 1);
    cw1 = cumsum(w1);
    wauc = sum(w1.*(cw(y==1)-cw1));
    sumw = cw1(length(cw1));
    sumw = sumw * (cw(length(cw)) - sumw);
    result = wauc / sumw;
end
end



    