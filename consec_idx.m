function [c_idx] = consec_idx(indices,threshold)

% Derives logical index of consecutive indices >= threshold
% e.g. if indices = [1 3 5 6 7 8 10 11 12 15], threshold = 3
%           c_idx = [0 0 1 1 1 1  1  1  1  0]
%     i.e. indices(c_idx) = [5 6 7 8 10 11 12]

k = [true,diff(indices)~=1];
s = cumsum(k);
x =  histc(s,1:s(end));
idx = find(k);
consecutive = idx(x>=threshold);

c_idx = false(size(indices));

for c=1:length(consecutive)
   x_idx = consecutive(c) == idx;
   c_idx(consecutive(c):(consecutive(c)+x(x_idx)-1)) = true;
end  