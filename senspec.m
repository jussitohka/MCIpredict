function [sen,spec,acc,bacc] = senspec(truelabels,predlabels,poslabel)

truelabels = truelabels == poslabel;
predlabels = predlabels == poslabel;

sen = length(find(predlabels == 1 & truelabels == 1))/length(find(truelabels == 1));
spec = length(find(predlabels == 0 & truelabels == 0))/length(find(truelabels == 0));
acc = length(find(predlabels == truelabels))/length(truelabels);
bacc = 0.5*sen + 0.5*spec;
end

