% [p_val, observeddifference, effectsize] = permTest_array(sample1, sample2, permutations [, varargin])
%
%       Permutation test (aka randomisation test), testing for a difference
%       in means between two samples. 
%
% In:
%       sample1 - array of sample 1 (e.g. peri-event signal; rows = trials)
%       sample2 - array of sample 2 (e.g. baseline signal; rows = trials)
%       permutations - the number of permutations
%
% Optional (name-value pairs):
%       sidedness - whether to test one- or two-sided:
%           'both' - test two-sided (default)
%           'smaller' - test one-sided, alternative hypothesis is that
%                       the mean of sample1 is smaller than the mean of
%                       sample2
%           'larger' - test one-sided, alternative hypothesis is that
%                      the mean of sample1 is larger than the mean of
%                      sample2
%       exact - whether or not to run an exact test, in which all possible
%               combinations are considered. this is only feasible for
%               relatively small sample sizes. the 'permutations' argument
%               will be ignored for an exact test. (1|0, default 0)
%       report_exact - report when exact test used (1|0, default: 1)
%       plotresult - whether or not to plot the distribution of randomised
%                    differences, along with the observed difference (1|0,
%                    default: 0)
%       showprogress - whether or not to show a progress bar. if 0, no bar
%                      is displayed; if showprogress > 0, the bar updates 
%                      every showprogress-th iteration.
%
% Out:  
%       p - the resulting p-value
%       observeddifference - the observed difference between the two
%                            samples, i.e. mean(sample1) - mean(sample2)
%       effectsize - the effect size
%
%  Copyright 2019 Philip Jean-Richard-dit-Bressel, UNSW Sydney
%  Based on Laurens R Krol 2019 permutationTest.m

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [p_val, observeddifference] = permTest_array(sample1, sample2, permutations, varargin)

% parsing input
p = inputParser;

addRequired(p, 'sample1', @isnumeric);
addRequired(p, 'sample2', @isnumeric);
addRequired(p, 'permutations', @isnumeric);

addParameter(p, 'sidedness', 'both', @(x) any(validatestring(x,{'both', 'smaller', 'larger'})));
addParameter(p, 'exact' , 0, @isnumeric);
addParameter(p, 'report_exact', 1, @isnumeric);
addParameter(p, 'showprogress', 0, @isnumeric);

parse(p, sample1, sample2, permutations, varargin{:})

sample1 = p.Results.sample1;
sample2 = p.Results.sample2;
permutations = p.Results.permutations;
sidedness = p.Results.sidedness;
exact = p.Results.exact;
report_exact = p.Results.report_exact;
showprogress = p.Results.showprogress;

allobservations = vertcat(sample1, sample2);
observeddifference = mean(sample1,1) - mean(sample2,1);
s1_n = size(sample1,1);
all_n = size(allobservations,1);
win_size = size(sample1,2);

w = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
if ~exact && permutations > nchoosek(all_n, s1_n)
    exact = 1;
    allcombinations = nchoosek(1:all_n, s1_n);
    if report_exact
      warning(['Number of permutations (%d) is higher than the number of possible combinations (%d);\n' ...
             'Running an exact test (minimum p = ' ...
             num2str(1/(nchoosek(all_n, s1_n)+1)) ')'], ...
             permutations, nchoosek(all_n, s1_n));
    end
    permutations = size(allcombinations, 1);
end
warning(w);

if showprogress, w = waitbar(0, 'Preparing test...', 'Name', 'permutationTest'); end

% running test
randomdifferences = zeros(win_size, permutations);
if showprogress, waitbar(0, w, sprintf('Permutation 1 of %d', permutations), 'Name', 'permutationTest'); end
for n = 1:permutations
    if showprogress && mod(n,showprogress) == 0, waitbar(n/permutations, w, sprintf('Permutation %d of %d', n, permutations)); end
    
    % selecting either next combination, or random permutation
    if exact, permutation = [allcombinations(n,:), setdiff(1:all_n, allcombinations(n,:))];
    else, permutation = randperm(all_n); end
    
    % dividing into two samples
    randomSample1 = allobservations(permutation(1:s1_n),:);
    randomSample2 = allobservations(permutation(s1_n+1:all_n),:);
    
    % saving differences between the two samples
    randomdifferences(:,n) = mean(randomSample1,1) - mean(randomSample2,1);
end
if showprogress, delete(w); end

% getting probability of finding observed difference from random permutations
p_val = zeros(1,win_size);

if strcmp(sidedness, 'both')
   for t = 1:win_size
   	p_val(t) = (sum(abs(randomdifferences(t,:)) > abs(observeddifference(t)))+1) / (permutations+1);
   end
elseif strcmp(sidedness, 'smaller')
    p = (length(find(randomdifferences < observeddifference))+1) / (permutations+1);
elseif strcmp(sidedness, 'larger')
    p = (length(find(randomdifferences > observeddifference))+1) / (permutations+1);
end

% % plotting result
% if plotresult
%     figure;
%     hist(randomdifferences, 20);
%     hold on;
%     xlabel('Random differences');
%     ylabel('Count')
%     od = plot(observeddifference, 0, '*r', 'DisplayName', sprintf('Observed difference.\nEffect size: %.2f,\np = %f', effectsize, p));
%     legend(od);
% end

end
