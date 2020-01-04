function bootCI = boot_CI(data,num_boots,sig)

%Bootstraps CI for mean population waveform

%%Input:
% data = data array (rows = trials, columns = time relative to event)
% num_boots = # bootstraps

%%Output:
% bootCI = LCI+UCI vector

%  Copyright 2019 Philip Jean-Richard-dit-Bressel, UNSW Sydney
%  Based on Colin Clifford 2018 bootstrap_CI.m 

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

[num_trials,window] = size(data);

% Minimum 2 trials (otherwise get funky signals due to inevitable crossing oscillations)
if num_trials > 2
   
   % Prep bootstrapping variables (one row for each bootstrap) ...
   data_boots = zeros(num_boots, window);
   bootCI = zeros(2,window);

   for b = 1:num_boots
      % bootstrap data + kernel across all trials ...
      trial_array = ceil((num_trials).*rand(1,num_trials));
      data_boots(b,:) = mean(data(trial_array,:));
   end
   
   %% Calculate bootstrap CI
   data_boots = sort(data_boots,1);

   lower_conf_index = ceil(num_boots*(sig/2))+1;
   upper_conf_index = floor(num_boots*(1-sig/2));

   bootCI(1,:) = data_boots(lower_conf_index,:);
   bootCI(2,:) = data_boots(upper_conf_index,:);

else
   fprintf('Less than 3 trials - bootstrapping skipped\n');
   bootCI = NaN;
end
