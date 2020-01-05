function [adjLCI,adjUCI] = CIadjust(LCI,UCI,mean_tr,n,adj_type)
   %% Adjusts CI according to type (adj_type)
      % Type = 1: extend each CI from mean by sqrt((n-1)/(n))
         % needs all inputs
      % Type = 2: expand CI by sqrt((n-1)/(n))
         % doesn't need mean
         
if adj_type == 1
   %% Extend CI from mean
   CI_fix = sqrt((n-1)/(n));
   fprintf(['CI extended from mean by ' num2str(1/CI_fix*100) 'pc\n']);
   adjUCI = (UCI-mean_tr)./CI_fix + mean_tr;
   adjLCI = mean_tr - (mean_tr-LCI)./CI_fix; 
   
elseif adj_type == 2
   %% Expand CI
   CI_fix = sqrt((n-1)/(n));
   fprintf(['CI expanded by ' num2str(1/CI_fix*100) 'pc\n']);
   CIchange = ((UCI - LCI)./CI_fix - (UCI - LCI))/2; %/2
   adjUCI = UCI+CIchange;
   adjLCI = LCI-CIchange;
end
end