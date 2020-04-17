%% ERT example
sig = .05;
consec_thresh = 340; % 1017.3Hz sample rate / 3Hz filter

% Graphing parameters
ylims = [-1 4];
xlims = [-2 5];
sig_level = linspace(3.5,4,7);
savefolder = ...
   'M:\Gavan McNally''s Lab\Staff folders\Philip';

%%
[n_Cp,ev_win] = size(ERT_test.Cp_off1);
[n_Cm,~] = size(ERT_test.Cm_off3);
timeline = linspace(-4,5,ev_win);

Cp_t_crit = tinv(1-sig/2,n_Cp-1);
Cm_t_crit = tinv(1-sig/2,n_Cm-1);

mean_Cp = mean(ERT_test.Cp_off1,1);
sem_Cp = sem(ERT_test.Cp_off1);
Cp_bCI = boot_CI(ERT_test.Cp_off1,1000,sig);
[adjLCI,adjUCI] = CIadjust(Cp_bCI(1,:),Cp_bCI(2,:),[],n_Cp,2);
Cp_bCIexp = [adjLCI;adjUCI];
Cp_tCI = [mean_Cp - sem_Cp*Cp_t_crit ; mean_Cp + sem_Cp*Cp_t_crit];

mean_Cm = mean(ERT_test.Cm_off3,1);
sem_Cm = sem(ERT_test.Cm_off3);
Cm_bCI = boot_CI(ERT_test.Cm_off3,1000,sig);
[adjLCI,adjUCI] = CIadjust(Cm_bCI(1,:),Cm_bCI(2,:),[],n_Cm,2);
Cm_bCIexp = [adjLCI;adjUCI];
Cm_tCI = [mean_Cm - sem_Cm*Cm_t_crit ; mean_Cm + sem_Cm*Cm_t_crit];

perm_p = permTest_array(ERT_test.Cp_off1,ERT_test.Cm_off3,1000);
diff_bCI = boot_diffCI(ERT_test.Cp_off1,ERT_test.Cm_off3,1000,sig);
[adjLCI,adjUCI] = CIadjust(diff_bCI(1,:),diff_bCI(2,:),[],n_Cm,2);
diff_bCIexp = [adjLCI;adjUCI];

%% Significance bars
Cp_bCIexp_sig = NaN(1,ev_win);
sig_idx = find((Cp_bCIexp(1,:) > 0) | (Cp_bCIexp(2,:) < 0));
consec = consec_idx(sig_idx,consec_thresh);
Cp_bCIexp_sig(sig_idx(consec)) = sig_level(3);

Cm_bCIexp_sig = NaN(1,ev_win);
sig_idx = find((Cm_bCIexp(1,:) > 0) | (Cm_bCIexp(2,:) < 0));
consec = consec_idx(sig_idx,consec_thresh);
Cm_bCIexp_sig(sig_idx(consec)) = sig_level(4);

Cp_tCI_sig = NaN(1,ev_win);
sig_idx = find((Cp_tCI(1,:) > 0) | (Cp_tCI(2,:) < 0));
consec = consec_idx(sig_idx,consec_thresh);
Cp_tCI_sig(sig_idx(consec)) = sig_level(5);

Cm_tCI_sig = NaN(1,ev_win);
sig_idx = find((Cm_tCI(1,:) > 0) | (Cm_tCI(2,:) < 0));
consec = consec_idx(sig_idx,consec_thresh);
Cm_tCI_sig(sig_idx(consec)) = sig_level(6);

perm_p_sig = NaN(1,ev_win);
sig_idx = find(perm_p < sig);
consec = consec_idx(sig_idx,consec_thresh);
perm_p_sig(sig_idx(consec)) = sig_level(1);

diff_bCIexp_sig = NaN(1,ev_win);
sig_idx = find((diff_bCIexp(1,:) > 0) | (diff_bCIexp(2,:) < 0));
consec = consec_idx(sig_idx,consec_thresh);
diff_bCIexp_sig(sig_idx(consec)) = sig_level(2);

diff_tCI_sig = NaN(1,ev_win);
sig_idx = ttest2(ERT_test.Cp_off1,ERT_test.Cm_off3);
sig_idx = find(sig_idx == 1);
consec = consec_idx(sig_idx,consec_thresh);
diff_tCI_sig(sig_idx(consec)) = sig_level(7);

%% Plot
figure; hold on
plot(timeline,mean_Cm,'Color',col_rep(3))
errorplot3(mean_Cm-sem_Cm,mean_Cm+sem_Cm,[-4 5],col_rep(3),.15)

plot(timeline,mean_Cp,'Color',col_rep(2))
errorplot3(mean_Cp-sem_Cp,mean_Cp+sem_Cp,[-4 5],col_rep(2),.15)

plot(timeline,Cp_bCIexp_sig,'Color',col_rep(2),'Marker','.')
text(xlims(1),sig_level(3),'CS+ bCI','Color',col_rep(2));
plot(timeline,Cm_bCIexp_sig,'Color',col_rep(3),'Marker','.')
text(xlims(1),sig_level(4),'CS- bCI','Color',col_rep(3));
plot(timeline,Cp_tCI_sig,'Color',col_rep(2),'Marker','.')
text(xlims(1),sig_level(5),'CS+ tCI','Color',col_rep(2));
plot(timeline,Cm_tCI_sig,'Color',col_rep(3),'Marker','.')
text(xlims(1),sig_level(6),'CS- tCI','Color',col_rep(3));
plot(timeline,perm_p_sig,'Color',col_rep(1),'Marker','.')
text(xlims(1),sig_level(1),'Perm','Color',col_rep(1));
plot(timeline,diff_bCIexp_sig,'Color',col_rep(4),'Marker','.')
text(xlims(1),sig_level(2),'Diff bCI','Color',col_rep(4));
plot(timeline,diff_tCI_sig,'Color',col_rep(4),'Marker','.')
text(xlims(1),sig_level(7),'Diff tCI','Color',col_rep(4));

plot([-0.5 -0.5],ylim,'k:')
plot(xlim,[0 0],'k--')

xlim(xlims);

print('-painters','-depsc2','-loose','-tiff', ...
	[savefolder '\ERTexampleResults.eps']);