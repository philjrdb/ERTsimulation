function SEM = sem(data)
% Standard error of mean - rows = trials

n = size(data,1);
SEM = std(data,0,1)/sqrt(n);