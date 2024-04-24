function [adjLCI,adjUCI] = CIadjust(LCI,UCI,est,n,adj_type)
%% Adjusts CI according to adj_type
   % Type 1 = extend CI from reference by sqrt((n)/(n-1))
      % needs all inputs
   % Type 2 = expand CI by sqrt((n)/(n-1))
      % doesn't need mean

%%Inputs:
% LCI = lower confidence limit vector
% UCI = upper confidence limit vector
% est = parameter estimate  vector (eg. mean, median) [only needed for adjustment type 1]
% n = experimental n
% adj_type = adjustment type (see above)

%%Output:
% adjLCI = adjusted LCI
% adjUCI = adjusted UCI

%% GNU
%  Copyright 2020 Philip Jean-Richard-dit-Bressel, UNSW Sydney

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

%%
CI_fix = sqrt((n)/(n-1));
if adj_type == 1
   %% Extend CI from reference
   fprintf(['CI extended from reference by ' num2str(CI_fix*100) 'pct\n']);
   adjUCI = (UCI-est).*CI_fix + est;
   adjLCI = est - (est-LCI).*CI_fix; 
   
elseif adj_type == 2
   %% Expand CI
   fprintf(['CI expanded by ' num2str(CI_fix*100) 'pct\n']);
   CIchange = ((UCI - LCI).*CI_fix - (UCI - LCI))/2; %/2
   adjUCI = UCI+CIchange;
   adjLCI = LCI-CIchange;
end
