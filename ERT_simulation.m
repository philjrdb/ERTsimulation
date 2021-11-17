%% ERT Simulation 
% Produces artificial signal (event-related transient; ERT) and control (null) populations datasets. 
% Uses Monte Carlo methods to assess Type I/II/III error rates for t confidence intervals (tCI), 
% bootstrapped confidence intervals (bCI) (n-factored expansion), and permutation tests.

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

%% Parameters
% Artificial dataset parameters 
% Peri-event window = vector of 1-100 (easy sum = %)
popn = 10000; % # signals (trials) per population
gen_subjs = 1000; % # subjects per population (0 = trial-based analysis)
max_trials = 30; % max # of trials/subject (random 1:max_trials+1)
sampling_rate = 10; % artificial sampling rate (# datapoints per 1 sec)
lowpass_filter = 2; % low-pass Hz to smooth dataset into waveforms
ERT = [0.1	0.6	0.95	1	0.8	0.5	0.3	0.2	0.1	0.05]; % ERT waveform (inserted at index 50 onwards)

% Analysis parameters
sample_sizes = [5, 8, 10, 15]; % sample sizes, eg. [5:10, 15, 20, 30, 40, 50, 100]
n_sim = 1000; % n Monte Carlo simulations per sample size
alpha = [5, 1]; % prescribed Type 1 error rate (% e.g. [5, 1])
consec_thresh = [3,5]; % consecutive threshold (only plots results of first 2)
resamples = 1000; % resamples for bootstrapping/permutation tests

% Misc
reload_boots = 1; %1 = try load bootstrap/permutation data from folder
savefolder = 'G:\Current\';

%% House-keeping
tic
folder = 'ERTsim';
mkdir(savefolder, folder);
savefolder = [savefolder folder];

a_leng = length(alpha);
max_n = sample_sizes(end);
consec_thresh_leng = length(consec_thresh);

% Idx fields in ERTsim structure (need to be reset between runs)
f_names = {'control_tCI_error' 'control_tCI_idx' 'control_true_tCI_error' 'control_true_tCI_idx' ...
           'signal_tCI_error' 'signal_tCI_idx' 'sig_CorRej_tCI' 'T3_tCI' ...       
           'control_bpCI_error' 'control_bpCI_idx' 'control_true_bpCI_error' 'control_true_bpCI_idx' ...
           'signal_bpCI_error' 'signal_bpCI_idx' 'sig_bpCI_CorRej' 'bpCI_T3' ...
           'control_bCIexp_error' 'control_bCIexp_idx' ...
           'control_true_bCIexp' 'control_true_bCIexp_idx' ...          
           'signal_bCIexp_error' 'signal_bCIexp_idx' 'fixsig_CorRej' 'fixT3'...
           'control_perm_idx' 'sig_CorRej_perm'};
f_leng = length(f_names);

%% Pre-allocate values/cells
for a = 1:a_leng
struct_name = ['alpha' num2str(alpha(a))];
sig = alpha(a)/100;

ERTsimsum.(struct_name).cont_bpCI_err{max_n} = [];
ERTsimsum.(struct_name).cont_bpCI_err_amount{max_n} = [];
ERTsimsum.(struct_name).cont_true_bpCI_err{max_n} = [];
ERTsimsum.(struct_name).cont_true_bpCI_err_amount{max_n} = [];
ERTsimsum.(struct_name).sig_bpCI_err{max_n} = [];
ERTsimsum.(struct_name).sig_CorRej_rate{max_n} = [];
ERTsimsum.(struct_name).T3{max_n} = [];

ERTsimsum.(struct_name).cont_bCIexp_err{max_n} = [];
ERTsimsum.(struct_name).cont_bCIexp_err_amount{max_n} = [];
ERTsimsum.(struct_name).cont_true_bCIexp_err{max_n} = [];
ERTsimsum.(struct_name).cont_true_bCIexp_err_amount{max_n} = [];
ERTsimsum.(struct_name).sig_bCIexp_err{max_n} = [];
ERTsimsum.(struct_name).fixsig_CorRej_rate{max_n} = [];
ERTsimsum.(struct_name).fixT3{max_n} = [];

for c = 1:consec_thresh_leng
   c_str = num2str(consec_thresh(c));
   
   ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err_amount']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err_amount']){max_n} = [];
   ERTsimsum.(struct_name).(['sig_bpCI_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['sig_CorRej_consec' c_str '_rate']){max_n} = [];
   ERTsimsum.(struct_name).(['T3_consec' c_str]){max_n} = [];

   ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err_amount']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err_amount']){max_n} = [];
   ERTsimsum.(struct_name).(['sig_bCIexp_consec' c_str '_err']){max_n} = [];
   ERTsimsum.(struct_name).(['fixsig_CorRej_consec' c_str '_rate']){max_n} = [];
   ERTsimsum.(struct_name).(['fixT3_consec' c_str]){max_n} = [];
end
end

%% Generate waveforms if control popn not in workspace   
if exist('control','var') ~= 1
   fprintf('\nGenerating population data\n');

   ERT_leng = length(ERT);

   control = awgn(zeros(popn,100),10);
   signal = zeros(popn,100);
   for s = 1:popn
      signal(s,50:50+ERT_leng-1) = ERT*abs(randn);
      if randn > 0
         rand_idx = randi(100-ERT_leng);
         control(s,rand_idx:(rand_idx+ERT_leng-1)) = ERT*abs(randn);
      end
   end

   signal = awgn(signal,10);

   lp_control = lowpass(control',lowpass_filter,sampling_rate)';
   lp_signal = lowpass(signal',lowpass_filter,sampling_rate)';

   mean_control = mean(lp_control);
   mean_signal = mean(lp_signal);

   close all
   figure; hold on
   plot(lp_control');
   plot(mean_control,'LineWidth',3);
   title(['Control population (n=' num2str(popn) ')'])
   saveas(figure(1), ...
      [savefolder '\Control population (n' num2str(popn) ').png']);
   saveas(figure(1), ...
      [savefolder '\Control population (n' num2str(popn) ').fig']);
   close all

   figure; hold on
   plot(lp_signal');
   plot(mean_signal,'LineWidth',3);
   title(['Signal population (n=' num2str(popn) ')']);
   saveas(figure(1), ...
      [savefolder '\Signal population (n' num2str(popn) ').png']);
   saveas(figure(1), ...
      [savefolder '\Signal population (n' num2str(popn) ').fig']);
   close all
else
   fprintf('\nUsing previously generated population data\n');
end

fprintf(['Control subject population: mean = ' num2str(mean(mean_control)) ', SD = ' num2str(std(mean_control)) ...
   '(range: ' num2str(min(mean_control)) ' - ' num2str(max(mean_control)) ')\n']);
fprintf(['Signal subject population: mean = ' num2str(mean(mean_signal)) ', SD = ' num2str(std(mean_signal)) ...
   '(range: ' num2str(min(mean_signal)) ' - ' num2str(max(mean_signal)) ')\n']);

% Subject population
if gen_subjs > 0
   if exist('subj_control','var') ~= 1
      fprintf('\nGenerating subject population data\n');
      subj_control = zeros(gen_subjs,100);
      subj_signal = zeros(gen_subjs,100);

      for n = 1:gen_subjs
         trials = ceil((popn).*rand(1,ceil(max_trials.*rand(1))+1));
         %fprintf(['Subj' num2str(n) ': ' num2str(length(trials)) ' trials\n']);
         subj_control(n,:) = mean(lp_control(trials,:));
         subj_signal(n,:) = mean(lp_signal(trials,:));
      end

      mean_control = mean(subj_control);
      mean_signal = mean(subj_signal);

      figure; hold on
      plot(subj_control');
      plot(mean_control,'LineWidth',3);
      title(['Control subject population (n=' num2str(gen_subjs) ')'])
      saveas(figure(1), ...
         [savefolder '\Control subject population (n' num2str(gen_subjs) ').png']);
      saveas(figure(1), ...
         [savefolder '\Control subject population (n' num2str(gen_subjs) ').fig']);
      close all

      figure; hold on
      plot(subj_signal');
      plot(mean_signal,'LineWidth',3);
      title(['Signal subject population (n=' num2str(gen_subjs) ')']);
      saveas(figure(1), ...
         [savefolder '\Signal subject population (n' num2str(gen_subjs) ').png']);
      saveas(figure(1), ...
         [savefolder '\Signal subject population (n' num2str(gen_subjs) ').fig']);
      close all

      fprintf(['Control population: mean = ' num2str(mean(mean_control)) ', SD = ' num2str(std(mean_control)) ...
         '(range: ' num2str(min(mean_control)) ' - ' num2str(max(mean_control)) ')\n']);
      fprintf(['Signal population: mean = ' num2str(mean(mean_signal)) ', SD = ' num2str(std(mean_signal)) ...
         '(range: ' num2str(min(mean_signal)) ' - ' num2str(max(mean_signal)) ')\n']);
   else
      fprintf('\nUsing previously generated subject population data\n');
   end
end

%% Loop through sample_sizes
sample_el = 1:length(sample_sizes);

for i_n = sample_el
   sample_n = sample_sizes(i_n);
   CI_fix = sqrt((sample_n-1)/(sample_n));
   n_struct = ['n' num2str(sample_n)];

   clc
   diary([savefolder '\ERTsim (n' ...
      num2str(sample_n) ').txt']);
   
   %% Load or bootstrap/permute
   if reload_boots == 1 && exist([savefolder '\ERTsim_' n_struct '.mat'],'file')
      fprintf('Loading previous Monte Carlo simulation\n'); 
      load([savefolder '\ERTsim_' n_struct '.mat']);      
   else
      clear ERTsim
      fprintf([' --- ERT Monte Carlo simulation ---\n' ...
            'Population = ' num2str(popn) ' waveforms per Control/Signal\n' ...
            num2str(gen_subjs) ' subjs with 1-' num2str(max_trials+1) ' waveforms\n']); 

      ERTsim.(n_struct)(n_sim) = struct('samples', [], 'perm',[], 'control_perm_p', [], 'signal_perm_p', [], ...
         'control_boots', zeros(resamples,100), 'signal_boots', zeros(resamples,100)); 
      for a = 1:a_leng
         struct_name = ['alpha' num2str(alpha(a))];
         alph = num2str(alpha(a));
      end
      
      for i = 1:n_sim
      if gen_subjs > 0 %subj-based
         ERTsim.(n_struct)(i).samples = ceil((gen_subjs).*rand(1,sample_n));
         ERTsim.(n_struct)(i).perm = ceil((gen_subjs).*rand(1,sample_n));
         
%          tic
%          fprintf('Bootstrap*resamples control+signal: ');
         ERTsim.(n_struct)(i).control_boots = boot_CI(subj_control(ERTsim.(n_struct)(i).samples,:),resamples,sig);
         ERTsim.(n_struct)(i).signal_boots = boot_CI(subj_signal(ERTsim.(n_struct)(i).samples,:),resamples,sig); 
%          toc
         
%          tic
%          fprintf('Permute*resamples control+signal: ');
         ERTsim.(n_struct)(i).control_perm_p = ...
            permTest_array(subj_control(ERTsim.(n_struct)(i).samples,:),subj_control(ERTsim.(n_struct)(i).perm,:),resamples,'report_exact',0);
         ERTsim.(n_struct)(i).signal_perm_p = ...
            permTest_array(subj_signal(ERTsim.(n_struct)(i).samples,:),subj_control(ERTsim.(n_struct)(i).perm,:),resamples,'report_exact',0);
%          toc         
         
      else % trial-based
         ERTsim.(n_struct)(i).samples = ceil((popn).*rand(1,sample_n));
         ERTsim.(n_struct)(i).perm = ceil((popn).*rand(1,sample_n));
         
         ERTsim.(n_struct)(i).control_boots = boot_CI(lp_control(ERTsim.(n_struct)(i).samples,:),resamples,sig);
         ERTsim.(n_struct)(i).signal_boots = boot_CI(lp_signal(ERTsim.(n_struct)(i).samples,:),resamples,sig); 

         ERTsim.(n_struct)(i).control_perm_p = ...
            permTest_array(subj_control(ERTsim.(n_struct)(i).samples,:),subj_control(ERTsim.(n_struct)(i).perm,:),resamples);
         ERTsim.(n_struct)(i).signal_perm_p = ...
            permTest_array(subj_signal(ERTsim.(n_struct)(i).samples,:),subj_control(ERTsim.(n_struct)(i).perm,:),resamples);
      end % subject vs trial
      end % i loop
   end % load or generate
   
   %% Assess results across time series      
   for a = 1:a_leng      
      struct_name = ['alpha' num2str(alpha(a))];
      alph = num2str(alpha(a)); 
      sig = alpha(a)/100;
      t_crit = tinv(1-sig/2,sample_n-1);
      
      fprintf(['\n - alpha = ' num2str(alpha(a)) ', n = ' num2str(sample_n) ' -\n' ...
         't-crit = ' num2str(t_crit) ', CI fix = divide range by ' num2str(CI_fix) ...
         ', ' num2str(consec_thresh) ' consecutive timepoints\n']); 

      for i = 1:n_sim
         mean_cont = mean(subj_control(ERTsim.(n_struct)(i).samples,:));
         sem_cont = std(subj_control(ERTsim.(n_struct)(i).samples,:))/sqrt(sample_n);
         mean_sig = mean(subj_signal(ERTsim.(n_struct)(i).samples,:));
         sem_sig = std(subj_signal(ERTsim.(n_struct)(i).samples,:))/sqrt(sample_n);    

         % Reset indices
         for f = 1:f_leng
            ERTsim.(n_struct)(i).([f_names{f} '_a' alph]) = []; 
         end
         
      for p = 1:100     
         ERTsim.(n_struct)(i).(['control_tUCI' '_a' alph])(p) = mean_cont(p) + t_crit*sem_cont(p);
         ERTsim.(n_struct)(i).(['control_tLCI' '_a' alph])(p) = mean_cont(p) - t_crit*sem_cont(p);
         ERTsim.(n_struct)(i).(['signal_tUCI' '_a' alph])(p) = mean_sig(p) + t_crit*sem_sig(p);
         ERTsim.(n_struct)(i).(['signal_tLCI' '_a' alph])(p) = mean_sig(p) - t_crit*sem_sig(p);
         
         ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p) = ...
            prctile(ERTsim.(n_struct)(i).control_boots(:,p),(100-(alpha(a)/2)));
         ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p) = ...
            prctile(ERTsim.(n_struct)(i).control_boots(:,p),(alpha(a)/2));
         ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p) = ...
            prctile(ERTsim.(n_struct)(i).signal_boots(:,p),(100-(alpha(a)/2)));
         ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p) = ...
            prctile(ERTsim.(n_struct)(i).signal_boots(:,p),(alpha(a)/2));

         %expand CI by undercoverage estimate
         CIchange = (((ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p) - ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p))/CI_fix)- ...
            (ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p) - ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p))); %/2
         ERTsim.(n_struct)(i).(['control_bUCIexp' '_a' alph])(p) = ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p)+CIchange;
         ERTsim.(n_struct)(i).(['control_bLCIexp' '_a' alph])(p) = ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p)-CIchange;

         CIchange = (((ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p) - ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p))/CI_fix)- ...
            (ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p) - ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p))); %/2
         ERTsim.(n_struct)(i).(['signal_bUCIexp' '_a' alph])(p) = ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p)+CIchange;
         ERTsim.(n_struct)(i).(['signal_bLCIexp' '_a' alph])(p) = ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p)-CIchange;
         
         %% Find control Type 1/coverage errors
         % permutation test
         if ERTsim.(n_struct)(i).control_perm_p(p) < sig
            ERTsim.(n_struct)(i).(['control_perm_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_perm_idx' '_a' alph]) p];
         end
            
         % tCI compared to popn mean
         if ERTsim.(n_struct)(i).(['control_tUCI' '_a' alph])(p) < mean_control(p)
            ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]) ...
               abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_tUCI' '_a' alph])(p))];
         end
         if ERTsim.(n_struct)(i).(['control_tLCI' '_a' alph])(p) > mean_control(p)
            ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]) ...
               abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_tLCI' '_a' alph])(p))];
         end
         
         % pCI compared to popn mean
         if ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p) < mean_control(p)
            ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph]) ...
               abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p))];
            if ERTsim.(n_struct)(i).(['control_bUCIexp' '_a' alph])(p) < mean_control(p)
               ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]) ...
                  abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_bUCIexp' '_a' alph])(p))];
            end
         end
         if ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p) > mean_control(p)
            ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph])...
               abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p))];
            if ERTsim.(n_struct)(i).(['control_bLCIexp' '_a' alph])(p) > mean_control(p)
               ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]) ...
                  abs(mean_control(p)-ERTsim.(n_struct)(i).(['control_bLCIexp' '_a' alph])(p))];
            end
         end
         
         % Compared to 0 (true null)
         if ERTsim.(n_struct)(i).(['control_tUCI' '_a' alph])(p) < 0
            ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]) ...
               abs(ERTsim.(n_struct)(i).(['control_tUCI' '_a' alph])(p))];
         end
         if ERTsim.(n_struct)(i).(['control_tLCI' '_a' alph])(p) > 0
            ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]) ...
               abs(ERTsim.(n_struct)(i).(['control_tLCI' '_a' alph])(p))];
         end
         if ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p) < 0
            ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]) ...
               abs(ERTsim.(n_struct)(i).(['control_bpUCI' '_a' alph])(p))];
            if ERTsim.(n_struct)(i).(['control_bUCIexp' '_a' alph])(p) < 0
               ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]) ...
                  abs(ERTsim.(n_struct)(i).(['control_bUCIexp' '_a' alph])(p))];
            end
         end
         if ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p) > 0
            ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]) ...
               abs(ERTsim.(n_struct)(i).(['control_bpLCI' '_a' alph])(p))];
            if (ERTsim.(n_struct)(i).(['control_bLCIexp' '_a' alph])(p) ) > 0
               ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]) = [ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]) ...
                  abs(ERTsim.(n_struct)(i).(['control_bLCIexp' '_a' alph])(p))];
            end
         end

         %% Signal coverage errors
         % tCI signal coverage errors
         if ERTsim.(n_struct)(i).(['signal_tUCI' '_a' alph])(p) < mean_signal(p)
            ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]) ...
               abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_tUCI' '_a' alph])(p))];
         end
         if ERTsim.(n_struct)(i).(['signal_tLCI' '_a' alph])(p) > mean_signal(p)
            ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]) ...
               abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_tLCI' '_a' alph])(p))];
         end         
         
         % pCI signal coverage errors
         if ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p) < mean_signal(p)
            ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]) ...
               abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p))];
            if ERTsim.(n_struct)(i).(['signal_bUCIexp' '_a' alph])(p) < mean_signal(p)
               ERTsim.(n_struct)(i).(['signal_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]) ...
                  abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_bUCIexp' '_a' alph])(p))];
            end
         end
         if ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p) > mean_signal(p)
            ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph]) p];
            ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]) ...
               abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p))];
            if ERTsim.(n_struct)(i).(['signal_bLCIexp' '_a' alph])(p) > mean_signal(p)
               ERTsim.(n_struct)(i).(['signal_bCIexp_idx' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bCIexp_idx' '_a' alph]) p];
               ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]) = [ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]) ...
                  abs(mean_signal(p)-ERTsim.(n_struct)(i).(['signal_bLCIexp' '_a' alph])(p))];
            end
         end
         if p > 49 && p < 50+ERT_leng
            if ERTsim.(n_struct)(i).signal_perm_p(p) < sig
               ERTsim.(n_struct)(i).(['sig_CorRej_perm' '_a' alph]) = [ERTsim.(n_struct)(i).(['sig_CorRej_perm' '_a' alph]) p];
            end
            if ERTsim.(n_struct)(i).(['signal_tLCI' '_a' alph])(p) > 0
               ERTsim.(n_struct)(i).(['sig_CorRej_tCI' '_a' alph]) = [ERTsim.(n_struct)(i).(['sig_CorRej_tCI' '_a' alph]) p];
            elseif ERTsim.(n_struct)(i).(['signal_tUCI' '_a' alph])(p) < 0
               ERTsim.(n_struct)(i).(['T3_tCI' '_a' alph]) = [ERTsim.(n_struct)(i).(['T3_tCI' '_a' alph]) p];
            end
            if ERTsim.(n_struct)(i).(['signal_bpLCI' '_a' alph])(p) > 0
               ERTsim.(n_struct)(i).(['sig_bpCI_CorRej' '_a' alph]) = [ERTsim.(n_struct)(i).(['sig_bpCI_CorRej' '_a' alph]) p];
            elseif ERTsim.(n_struct)(i).(['signal_bpUCI' '_a' alph])(p) < 0
               ERTsim.(n_struct)(i).(['bpCI_T3' '_a' alph]) = [ERTsim.(n_struct)(i).(['bpCI_T3' '_a' alph]) p];
            end
            if ERTsim.(n_struct)(i).(['signal_bLCIexp' '_a' alph])(p) > 0
               ERTsim.(n_struct)(i).(['fixsig_CorRej' '_a' alph]) = [ERTsim.(n_struct)(i).(['fixsig_CorRej' '_a' alph]) p];
            elseif ERTsim.(n_struct)(i).(['signal_bUCIexp' '_a' alph])(p) < 0
               ERTsim.(n_struct)(i).(['fixT3' '_a' alph]) = [ERTsim.(n_struct)(i).(['fixT3' '_a' alph]) p];
            end
         end
      end

      %% Isolate consecutively significant timepoints
      for c = 1:consec_thresh_leng
         c_str = num2str(consec_thresh(c));
         
         % Permutation
         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_perm_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_perm_p_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_perm_idx' '_a' alph])(c_idx);
         c_idx = consec_idx(ERTsim.(n_struct)(i).(['sig_CorRej_perm' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['sig_CorRej_perm_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['sig_CorRej_perm' '_a' alph])(c_idx);

         % Control 
         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_tCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_tCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_tCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_bpCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bpCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_bpCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_bCIexp_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_bCIexp_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph])(c_idx);

         % True null
         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_true_tCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_tCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_true_tCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_true_bpCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_bpCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_true_bpCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['control_true_bCIexp_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_bCIexp_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['control_true_bCIexp_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph])(c_idx);

         % Signal
         c_idx = consec_idx(ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['signal_tCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['signal_tCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['signal_tCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['signal_bpCI_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['signal_bpCI_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['signal_bpCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['signal_bCIexp_idx_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bCIexp_idx' '_a' alph])(c_idx);
         ERTsim.(n_struct)(i).(['signal_bCIexp_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['sig_CorRej_tCI' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['sig_CorRej_tCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['sig_CorRej_tCI' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['sig_bpCI_CorRej' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['sig_CorRej_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['sig_bpCI_CorRej' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['fixsig_CorRej' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['fixsig_CorRej_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['fixsig_CorRej' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['T3_tCI' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['T3_tCI_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['T3_tCI' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['bpCI_T3' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['T3_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['bpCI_T3' '_a' alph])(c_idx);

         c_idx = consec_idx(ERTsim.(n_struct)(i).(['fixT3' '_a' alph]),consec_thresh(c));
         ERTsim.(n_struct)(i).(['fixT3_consec' c_str '_a' alph]) = ERTsim.(n_struct)(i).(['fixT3' '_a' alph])(c_idx);
      end

      %% Error/win rates
      ERTsim.(n_struct)(i).(['control_perm_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_perm_idx' '_a' alph]));
      ERTsim.(n_struct)(i).(['sig_CorRej_perm_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_CorRej_perm' '_a' alph]))/ERT_leng*100;
      
      ERTsim.(n_struct)(i).(['control_tCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_tCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_tCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_tCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_tCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_tCI_error' '_a' alph]));
      
      ERTsim.(n_struct)(i).(['control_bpCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_bpCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_bpCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_bpCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_bpCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_bpCI_error' '_a' alph]));

      ERTsim.(n_struct)(i).(['control_bCIexp_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_bCIexp_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_bCIexp_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_bCIexp_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]));
      ERTsim.(n_struct)(i).(['control_true_bCIexp_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_bCIexp' '_a' alph]));
      
      ERTsim.(n_struct)(i).(['signal_tCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['signal_tCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_tCI_error' '_a' alph]));

      ERTsim.(n_struct)(i).(['signal_bpCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['signal_bpCI_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_bpCI_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['signal_bCIexp_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]));
      ERTsim.(n_struct)(i).(['signal_bCIexp_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_bCIexp_error' '_a' alph]));

      ERTsim.(n_struct)(i).(['sig_CorRej_tCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_CorRej_tCI' '_a' alph]))/ERT_leng*100;
      ERTsim.(n_struct)(i).(['T3_tCI_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['T3_tCI' '_a' alph]))/ERT_leng*100;

      ERTsim.(n_struct)(i).(['sig_CorRej_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_bpCI_CorRej' '_a' alph]))/ERT_leng*100;
      ERTsim.(n_struct)(i).(['T3_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['bpCI_T3' '_a' alph]))/ERT_leng*100;
      ERTsim.(n_struct)(i).(['fixsig_CorRej_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['fixsig_CorRej' '_a' alph]))/ERT_leng*100;
      ERTsim.(n_struct)(i).(['fixT3_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['fixT3' '_a' alph]))/ERT_leng*100;
      
      for c = 1:consec_thresh_leng
         c_str = num2str(consec_thresh(c));
         
         ERTsim.(n_struct)(i).(['control_perm_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_perm_p_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['sig_CorRej_perm_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_CorRej_perm_consec' c_str '_a' alph]))/ERT_leng*100;

         ERTsim.(n_struct)(i).(['signal_tCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_tCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['signal_tCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_tCI_consec' c_str '_a' alph]));

         ERTsim.(n_struct)(i).(['control_tCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_tCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_tCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_tCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_tCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_tCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_tCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_tCI_consec' c_str '_a' alph]));

         ERTsim.(n_struct)(i).(['control_bpCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_bpCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_bpCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_bpCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_bpCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_bpCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_bpCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_bpCI_consec' c_str '_a' alph]));

         ERTsim.(n_struct)(i).(['control_bCIexp_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_bCIexp_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_bCIexp_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_bCIexp_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_bCIexp_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['control_true_bCIexp_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['control_true_bCIexp_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['control_true_bCIexp_consec' c_str '_a' alph]));

         ERTsim.(n_struct)(i).(['signal_bpCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_bpCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['signal_bpCI_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_bpCI_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['signal_bCIexp_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['signal_bCIexp_consec' c_str '_a' alph]));
         ERTsim.(n_struct)(i).(['signal_bCIexp_consec' c_str '_amount' '_a' alph]) = mean(ERTsim.(n_struct)(i).(['signal_bCIexp_consec' c_str '_a' alph]));

         ERTsim.(n_struct)(i).(['sig_CorRej_tCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_CorRej_tCI_consec' c_str '_a' alph]))/ERT_leng*100;
         ERTsim.(n_struct)(i).(['T3_tCI_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['T3_tCI_consec' c_str '_a' alph]))/ERT_leng*100;

         ERTsim.(n_struct)(i).(['sig_CorRej_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['sig_CorRej_consec' c_str '_a' alph]))/ERT_leng*100;
         ERTsim.(n_struct)(i).(['T3_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['T3_consec' c_str '_a' alph]))/ERT_leng*100;
         ERTsim.(n_struct)(i).(['fixsig_CorRej_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['fixsig_CorRej_consec' c_str '_a' alph]))/ERT_leng*100;
         ERTsim.(n_struct)(i).(['fixT3_consec' c_str '_rate' '_a' alph]) = length(ERTsim.(n_struct)(i).(['fixT3_consec' c_str '_a' alph]))/ERT_leng*100;
      end %
      end % end i loop

      %% Aggregate results
      ERTsimsum.(struct_name).control_perm_rate{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_perm_rate' '_a' alph]));
      ERTsimsum.(struct_name).sig_CorRej_perm_rate{sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_perm_rate' '_a' alph]));
      ERTsimsum.(struct_name).control_perm_FWER{sample_n} = (ERTsimsum.(struct_name).control_perm_rate{sample_n} > 0);
 
      ERTsimsum.(struct_name).cont_tCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_tCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_tCI_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_tCI_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_tCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_tCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_tCI_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_tCI_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_tCI_FWER{sample_n} = (ERTsimsum.(struct_name).cont_tCI_err{sample_n} > 0);
      ERTsimsum.(struct_name).cont_true_tCI_FWER{sample_n} = (ERTsimsum.(struct_name).cont_true_tCI_err{sample_n} > 0);

      ERTsimsum.(struct_name).sig_tCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_tCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).sig_CorRej_tCI_rate{sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_tCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).T3_t{sample_n} = vertcat(ERTsim.(n_struct)(:).(['T3_tCI_rate' '_a' alph]));

      ERTsimsum.(struct_name).cont_bpCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bpCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_bpCI_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bpCI_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_bpCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bpCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_bpCI_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bpCI_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_bpCI_FWER{sample_n} = (ERTsimsum.(struct_name).cont_bpCI_err{sample_n} > 0);
      ERTsimsum.(struct_name).cont_true_bpCI_FWER{sample_n} = (ERTsimsum.(struct_name).cont_true_bpCI_err{sample_n} > 0);
 
      ERTsimsum.(struct_name).sig_bpCI_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_bpCI_rate' '_a' alph]));
      ERTsimsum.(struct_name).sig_CorRej_rate{sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_rate' '_a' alph]));
      ERTsimsum.(struct_name).T3{sample_n} = vertcat(ERTsim.(n_struct)(:).(['T3_rate' '_a' alph]));

      ERTsimsum.(struct_name).cont_bCIexp_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bCIexp_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_bCIexp_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bCIexp_rate' '_a' alph]));
      ERTsimsum.(struct_name).cont_bCIexp_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bCIexp_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_true_bCIexp_err_amount{sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bCIexp_amount' '_a' alph]));
      ERTsimsum.(struct_name).cont_bCIexp_FWER{sample_n} = (ERTsimsum.(struct_name).cont_bCIexp_err{sample_n} > 0);
      ERTsimsum.(struct_name).cont_true_bCIexp_FWER{sample_n} = (ERTsimsum.(struct_name).cont_true_bCIexp_err{sample_n} > 0);

      ERTsimsum.(struct_name).sig_bCIexp_err{sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_bCIexp_rate' '_a' alph]));
      ERTsimsum.(struct_name).fixsig_CorRej_rate{sample_n} = vertcat(ERTsim.(n_struct)(:).(['fixsig_CorRej_rate' '_a' alph]));
      ERTsimsum.(struct_name).fixT3{sample_n} = vertcat(ERTsim.(n_struct)(:).(['fixT3_rate' '_a' alph]));
           
      for c = 1:consec_thresh_leng
         c_str = num2str(consec_thresh(c));
         
         ERTsimsum.(struct_name).(['control_perm_consec' c_str '_rate']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_perm_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' c_str '_rate']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_perm_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['control_perm_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['control_perm_consec' c_str '_rate']){sample_n} > 0);

         ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_tCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_tCI_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_tCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_tCI_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_err']){sample_n} > 0);
         ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_err']){sample_n} > 0);

         ERTsimsum.(struct_name).(['sig_tCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_tCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' c_str '_rate']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_tCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['T3_tCI_consec' c_str]){sample_n} = vertcat(ERTsim.(n_struct)(:).(['T3_tCI_consec' c_str '_rate' '_a' alph]));

         ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bpCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bpCI_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bpCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bpCI_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err']){sample_n} > 0);
         ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err']){sample_n} > 0);

         ERTsimsum.(struct_name).(['sig_bpCI_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_bpCI_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['sig_CorRej_consec' c_str '_rate']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['sig_CorRej_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['T3_consec' c_str]){sample_n} = vertcat(ERTsim.(n_struct)(:).(['T3_consec' c_str '_rate' '_a' alph]));

         ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bCIexp_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bCIexp_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_bCIexp_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err_amount']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['control_true_bCIexp_consec' c_str '_amount' '_a' alph]));
         ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err']){sample_n} > 0);
         ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_FWER']){sample_n} = (ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err']){sample_n} > 0);

         ERTsimsum.(struct_name).(['sig_bCIexp_consec' c_str '_err']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['signal_bCIexp_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['fixsig_CorRej_consec' c_str '_rate']){sample_n} = vertcat(ERTsim.(n_struct)(:).(['fixsig_CorRej_consec' c_str '_rate' '_a' alph]));
         ERTsimsum.(struct_name).(['fixT3_consec' c_str]){sample_n} = vertcat(ERTsim.(n_struct)(:).(['fixT3_consec' c_str '_rate' '_a' alph]));
      end
      
      %% Report results
      fprintf('\n - Permutation test results -\n');
      fprintf(['Avg Type I rate: ' num2str(mean(ERTsimsum.(struct_name).control_perm_rate{sample_n})) ' (' ...
         num2str(min(ERTsimsum.(struct_name).control_perm_rate{sample_n})) ...
         '-' num2str(max(ERTsimsum.(struct_name).control_perm_rate{sample_n})) ')\n']);
      fprintf(['FWER: ' num2str(mean(ERTsimsum.(struct_name).control_perm_FWER{sample_n})*100) 'pc\n']);
      fprintf(['Avg signal correct reject rate: ' num2str(mean(ERTsimsum.(struct_name).sig_CorRej_perm_rate{sample_n})) ' (' ...
         num2str(min(ERTsimsum.(struct_name).sig_CorRej_perm_rate{sample_n})) ...
         '-' num2str(max(ERTsimsum.(struct_name).sig_CorRej_perm_rate{sample_n})) ')\n']);

      fprintf('\n - Permutation consec results -\n');
      fprintf(['Avg Type I rate: ' num2str(mean(ERTsimsum.(struct_name).(['control_perm_consec' c_str '_rate']){sample_n})) ' (' ...
         num2str(min(ERTsimsum.(struct_name).(['control_perm_consec' c_str '_rate']){sample_n})) ...
         '-' num2str(max(ERTsimsum.(struct_name).(['control_perm_consec' c_str '_rate']){sample_n})) ')\n']);
      fprintf(['FWER: ' num2str(mean(ERTsimsum.(struct_name).(['control_perm_consec' c_str '_FWER']){sample_n})*100) 'pc\n']);
      fprintf(['Avg signal correct reject rate: ' num2str(mean(ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' c_str '_rate']){sample_n})) ' (' ...
         num2str(min(ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' c_str '_rate']){sample_n})) ...
         '-' num2str(max(ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' c_str '_rate']){sample_n})) ')\n']);
      
      fprintf('\n - Standard-t coverage error rates -\n');
      ERT_sim_report(ERTsimsum.(struct_name).cont_tCI_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_tCI_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_tCI_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_tCI_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).sig_tCI_err{sample_n}, ...
         ERTsimsum.(struct_name).sig_CorRej_tCI_rate{sample_n}, ...
         ERTsimsum.(struct_name).T3_t{sample_n});

      fprintf('\n - bpCI coverage error rates -\n');
      ERT_sim_report(ERTsimsum.(struct_name).cont_bpCI_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_bpCI_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_bpCI_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_bpCI_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).sig_bpCI_err{sample_n}, ...
         ERTsimsum.(struct_name).sig_CorRej_rate{sample_n}, ...
         ERTsimsum.(struct_name).T3{sample_n});
      
      fprintf('\n - bCIexp coverage error rates -\n');
      ERT_sim_report(ERTsimsum.(struct_name).cont_bCIexp_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_bCIexp_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_bCIexp_err{sample_n}, ...
         ERTsimsum.(struct_name).cont_true_bCIexp_err_amount{sample_n}, ...
         ERTsimsum.(struct_name).sig_bCIexp_err{sample_n}, ...
         ERTsimsum.(struct_name).fixsig_CorRej_rate{sample_n}, ...
         ERTsimsum.(struct_name).fixT3{sample_n});
      
      for c = 1:consec_thresh_leng
         c_str = num2str(consec_thresh(c));
         
         fprintf(['\n - Consec' c_str ' t coverage error rates -\n']);
         ERT_sim_report(ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_tCI_consec' c_str '_err_amount']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_true_tCI_consec' c_str '_err_amount']){sample_n}, ...
            ERTsimsum.(struct_name).(['sig_tCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' c_str '_rate']){sample_n}, ...
            ERTsimsum.(struct_name).(['T3_tCI_consec' c_str]){sample_n});

         fprintf(['\n - bpCI consec' c_str ' coverage error rates -\n']);
         ERT_sim_report(ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_bpCI_consec' c_str '_err_amount']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_true_bpCI_consec' c_str '_err_amount']){sample_n}, ...
            ERTsimsum.(struct_name).(['sig_bpCI_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['sig_CorRej_consec' c_str '_rate']){sample_n}, ...
            ERTsimsum.(struct_name).(['T3_consec' c_str]){sample_n});

         fprintf(['\n - bCIexp consec' c_str ' coverage error rates -\n']);
         ERT_sim_report(ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_bCIexp_consec' c_str '_err_amount']){sample_n},...
            ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['cont_true_bCIexp_consec' c_str '_err_amount']){sample_n}, ...
            ERTsimsum.(struct_name).(['sig_bCIexp_consec' c_str '_err']){sample_n}, ...
            ERTsimsum.(struct_name).(['fixsig_CorRej_consec' c_str '_rate']){sample_n}, ...
            ERTsimsum.(struct_name).(['fixT3_consec' c_str]){sample_n}); 
      end
      
      diary off
   end
   
   save([savefolder '\ERTsim_' n_struct '.mat'], 'ERTsim');
   clear ERTsim;
   
end

%% Plot non-coverage rates
for a = 1:a_leng
   struct_name = ['alpha' num2str(alpha(a))];
   sig = alpha(a)/100;

   figure; hold on
   t1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).cont_tCI_err(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).cont_tCI_err(sample_sizes)), ...
      'k:x', 'LineWidth', 1.5);
   nom = plot(xlim,[alpha(a) alpha(a)],'k--');
   t2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(1)) '_err'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(1)) '_err'])(sample_sizes)), ...
      'k--d', 'LineWidth', 1.5);
   t3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(2)) '_err'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(2)) '_err'])(sample_sizes)), ...
      'k-o', 'LineWidth', 1.5);
   b1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).cont_bCIexp_err(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).cont_bCIexp_err(sample_sizes)), ...
      'r:x', 'LineWidth', 1.5);
   b2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(1)) '_err'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(1)) '_err'])(sample_sizes)), ...
      'r--d', 'LineWidth', 1.5);   
   b3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(2)) '_err'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(2)) '_err'])(sample_sizes)), ...
      'r-o', 'LineWidth', 1.5);
   perm1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).control_perm_rate(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).control_perm_rate(sample_sizes)), ...
      'b:x', 'LineWidth', 1.5);
   perm2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'b--d', 'LineWidth', 1.5);
   perm3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'b-o', 'LineWidth', 1.5);

   % Legends etc
   legend([nom t1 t2 t3 b1 b2 b3 perm1 perm2 perm3],...
      {'Nominal rate', ...
      'Standard tCI', ['Consec' num2str(consec_thresh(1)) ' tCI'], ['Consec' num2str(consec_thresh(2)) ' tCI'],...
      'bCIexp', ['Consec' num2str(consec_thresh(1)) ' bCIexp'], ['Adj Consec' num2str(consec_thresh(2)) ' bCIexp'], ...
      'Perm Type I', ['Perm consec'  num2str(consec_thresh(1))], ['Perm consec'  num2str(consec_thresh(2))]}, ...
      'Location', 'best');
   ylabel('Non-coverage error rate (%)','interpreter','none')
   ylims = ylim;
   ylims(1) = 0;
   ylim(ylims);
   xlabel('Sample n');
   xticks(sample_el);
   xticklabels(string(sample_sizes));
   title(['Non-coverage error rate (a = ' num2str(sig) ')']);
   set(gcf,'Position',get(0,'Screensize'));

   saveas(figure(1),[savefolder '\CI non-coverage rate (' struct_name ').png']);
   saveas(figure(1),[savefolder '\CI non-coverage rate (' struct_name ').fig']);
   close all
   
   %% FWER
   figure; hold on
   t1 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).cont_tCI_FWER(sample_sizes))*100, ...
      'k:x', 'LineWidth', 1.5);
   nom = plot(xlim,[alpha(a) alpha(a)],'k--');
   t2 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(1)) '_FWER'])(sample_sizes))*100, ...
      'k--d', 'LineWidth', 1.5);
   t3 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_tCI_consec' num2str(consec_thresh(2)) '_FWER'])(sample_sizes))*100, ...
      'k-o', 'LineWidth', 1.5);
   b1 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).cont_bCIexp_FWER(sample_sizes))*100, ...
      'r:x', 'LineWidth', 1.5);
   b2 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(1)) '_FWER'])(sample_sizes))*100, ...
      'r--d', 'LineWidth', 1.5);
   b3 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['cont_bCIexp_consec' num2str(consec_thresh(2)) '_FWER'])(sample_sizes))*100, ...
      'r-o', 'LineWidth', 1.5);
   perm1 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).control_perm_FWER(sample_sizes))*100, ...
      'b:x', 'LineWidth', 1.5);
   perm2 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(1)) '_FWER'])(sample_sizes))*100, ...
      'b--d', 'LineWidth', 1.5);
   perm3 = plot(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['control_perm_consec' num2str(consec_thresh(2)) '_FWER'])(sample_sizes))*100, ...
      'b-o', 'LineWidth', 1.5);

   % Legends etc
   legend([nom t1 t2 t3 b1 b2 b3 perm1 perm2 perm3],...
      {'Nominal rate', ...
      'Standard tCI', ['Consec' num2str(consec_thresh(1)) ' tCI'], ['Consec' num2str(consec_thresh(2)) ' tCI'],...
      'bCIexp', ['Consec' num2str(consec_thresh(1)) ' bCIexp'], ['Consec' num2str(consec_thresh(2)) ' bCIexp'], ...
      'Perm Type I', ['Perm consec'  num2str(consec_thresh(1))], ['Perm consec'  num2str(consec_thresh(2))]}, ...
      'Location', 'best');
   ylabel('FWER (%)','interpreter','none')
   ylims = ylim;
   ylims(1) = 0;
   ylim(ylims);
   xlabel('Sample n');
   xticks(sample_el);
   xticklabels(string(sample_sizes));
   title(['FWER (a = ' num2str(sig) ')']);
   set(gcf,'Position',get(0,'Screensize'));

   saveas(figure(1),[savefolder '\FWER (' struct_name ').png']);
   saveas(figure(1),[savefolder '\FWER (' struct_name ').fig']);
   close all

   %% Plot Type II rates
   figure; hold on
   t1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).sig_CorRej_tCI_rate(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).sig_CorRej_tCI_rate(sample_sizes)), ...
      'k:x', 'LineWidth', 1.5);
   t2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'k--d', 'LineWidth', 1.5);
   t3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'k-o', 'LineWidth', 1.5);
   b1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).fixsig_CorRej_rate(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).fixsig_CorRej_rate(sample_sizes)), ...
      'r:x', 'LineWidth', 1.5);
   b2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'r--d', 'LineWidth', 1.5);
   b3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'r-o', 'LineWidth', 1.5);
   perm1 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).sig_CorRej_perm_rate(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).sig_CorRej_perm_rate(sample_sizes)), ...
      'b:x', 'LineWidth', 1.5);
   perm2 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'b--d', 'LineWidth', 1.5);
   perm3 = errorbar(sample_el, ...
      cellfun(@mean,ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      cellfun(@(x) std(x)/sqrt(n_sim),ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'b-o', 'LineWidth', 1.5);

   % Legends etc
   legend([t1 t2 t3 b1 b2 b3 perm1 perm2 perm3],...
      {'Standard tCI', ['Consec' num2str(consec_thresh(1)) ' tCI'], ['Consec' num2str(consec_thresh(2)) ' tCI'],...
      'bCIexp', ['Consec' num2str(consec_thresh(1)) ' bCIexp'], ['Consec' num2str(consec_thresh(2)) ' bCIexp'], ...
      'Perm Type I', ['Perm consec'  num2str(consec_thresh(1))], ['Perm consec'  num2str(consec_thresh(2))]}, ...
      'Location', 'best');
   ylabel('Correct reject rate (%)','interpreter','none')
   xlabel('Sample n');
   xticks(sample_el);
   xticklabels(string(sample_sizes));
   title(['Correct reject rate (a = ' num2str(sig) ')']);
   set(gcf,'Position',get(0,'Screensize'));

   saveas(figure(1),[savefolder '\Correct reject rate (' struct_name ').png']);
   saveas(figure(1),[savefolder '\Correct reject rate (' struct_name ').fig']);
   close all

   fprintf([struct_name ' done. ']);
   toc
   
   %% Plot FW Type II rate
   figure; hold on
   t1 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).sig_CorRej_tCI_rate(sample_sizes)), ...
      'k:x', 'LineWidth', 1.5);
   t2 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'k--d', 'LineWidth', 1.5);
   t3 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['sig_CorRej_tCI_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'k-o', 'LineWidth', 1.5);
   b1 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).fixsig_CorRej_rate(sample_sizes)), ...
      'r:x', 'LineWidth', 1.5);
   b2 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'r--d', 'LineWidth', 1.5);
   b3 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['fixsig_CorRej_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'r-o', 'LineWidth', 1.5);
   perm1 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).sig_CorRej_perm_rate(sample_sizes)), ...
      'b:x', 'LineWidth', 1.5);
   perm2 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(1)) '_rate'])(sample_sizes)), ...
      'b--d', 'LineWidth', 1.5);
   perm3 = plot(sample_el, ...
      cellfun(@(x) (100-(sum(x > 0)/n_sim*100)),ERTsimsum.(struct_name).(['sig_CorRej_perm_consec' num2str(consec_thresh(2)) '_rate'])(sample_sizes)), ...
      'b-o', 'LineWidth', 1.5);

   % Legends etc
   legend([t1 t2 t3 b1 b2 b3 perm1 perm2 perm3],...
      {'tCI', ['Consec' num2str(consec_thresh(1)) ' tCI'], ['Consec' num2str(consec_thresh(2)) ' tCI'],...
      'bCIexp', ['Consec' num2str(consec_thresh(1)) ' bCIexp'], ['Adj Consec' num2str(consec_thresh(2)) ' bCIexp'], ...
      'Perm Type I', ['Perm consec'  num2str(consec_thresh(1))], ['Perm consec'  num2str(consec_thresh(2))]}, ...
      'Location', 'best');
   ylabel('FW Type II rate (%)','interpreter','none')
   xlabel('Sample n');
   xticks(sample_el);
   xticklabels(string(sample_sizes));
   title(['FW Type II (a = ' num2str(sig) ')']);
   set(gcf,'Position',get(0,'Screensize'));

   saveas(figure(1),[savefolder '\FW Type II (' struct_name ').png']);
   saveas(figure(1),[savefolder '\FW Type II (' struct_name ').fig']);
   close all

   fprintf([struct_name ' done. ']);
   toc
end

save([savefolder '\' folder '.mat'], ...
   'ERTsimsum', 'control', 'signal', 'lp_control', 'lp_signal', 'subj_control', 'subj_signal');

fprintf([folder ' done + saved. ']);
toc
