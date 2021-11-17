function ERT_sim_report(cont_err,cont_err_amount,cont_true_err,cont_true_err_amount, ...
   sig_err,sig_win_rate,T3)
   
% Reports results of simulation for null/control (cont) and ERT/signal (sig))

%% By Philip Jean-Richard-dit-Bressel, UNSW Sydney, 2020
% Feel free to use with citation: Jean-Richard-dit-Bressel et al. (2020). https://doi.org/10.3389/fnmol.2020.00014

%% GNU
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

fprintf(['Avg control error rate: ' num2str(mean(cont_err)) ' (' ...
   num2str(min(cont_err)) '-' num2str(max(cont_err)) ')\n']);
fprintf(['Control FWER: ' num2str(mean(cont_err > 0)*100) 'pc\n']);
fprintf(['Avg control error amount: ' num2str(mean(cont_err_amount,'omitnan')) ' (' ...
   num2str(min(cont_err_amount)) '-' num2str(max(cont_err_amount)) ')\n']);
fprintf(['Avg true error rate: ' num2str(mean(cont_true_err)) ' (' ...
   num2str(min(cont_true_err)) '-' num2str(max(cont_true_err)) ')\n']);
fprintf(['True control FWER: ' num2str(mean(cont_true_err > 0)*100) 'pc\n']);
fprintf(['Avg true error amount: ' num2str(mean(cont_true_err_amount,'omitnan')) ' (' ...
   num2str(min(cont_true_err_amount)) '-' num2str(max(cont_true_err_amount)) ')\n']);
fprintf(['Avg signal error rate: ' num2str(mean(sig_err)) ' (' ...
   num2str(min(sig_err)) '-' num2str(max(sig_err)) ')\n']);
fprintf(['Signal FWER: ' num2str(mean(sig_err > 0)*100) 'pc\n']);
fprintf(['Avg signal correct reject rate: ' num2str(mean(sig_win_rate)) ' (' ...
   num2str(min(sig_win_rate)) '-' num2str(max(sig_win_rate)) ')\n']);
fprintf(['FW correct reject rate: ' num2str(mean(sig_win_rate > 0)*100) 'pc\n']);
fprintf(['Avg signal T3 error rate: ' num2str(mean(T3)) ' (' ...
   num2str(min(T3)) '-' num2str(max(T3)) ')\n']);
