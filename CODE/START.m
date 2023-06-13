% Joint-acceleration-model for association of multiple coordinates, with
% maximum likelihood and MCMC estimators, outlier detection and adjustment,
% Kalman filter, Kalman smoother, particle filter, particle smoother, and
% MAP smoother (circa 1500 lines of code all together)

close all
clear all
clc

%%%%%%%%%%%%%%%
del = 1/15*4;   % - every 5th frame (refreshing 15 frame per second)
%%%%%%%%%%%%%%%

%%%%%%%%%%%%
nop = 100;   % - number of particles
%%%%%%%%%%%%

%%%%%%%%%%%%%%
resx = 1280;   % - pixel resolution along x coordinate
resy = 720;    % - pixel resolution along y coordinate
%%%%%%%%%%%%%% %   (these inputs are very important for evaluation of
               %    likelihood and computing weights in particle filter) 

%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'output_1187'; % or 1187 or 1173 or 95 or 1126
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs = struct;%%%%%%%% insert 0 (NO), 1 (YES) or [] (DECIDE ONLINE) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.plot_raw = 0;     % - 1 to plot raw coordinates
prefs.plot_trun = 1;    % - 1 to plot truncated coordinates
%------------------------
prefs.adj_dums = 1;     % - 1 to adjust outliers using dummies and mle
prefs.adj_choice = 1;   % - 1 to adjust using dummies, 0 using means, 
                        %   [] to decide online after viewing adjustments
                        %   (it is better to choose means for longer time
                        %   series since estimating model with dummies may
                        %   require considerable time before convergence
%------------------------
prefs.mcmc_ests = 1;    % - 1 to perform MCMC estimation
prefs.plot_dgns = 0;    % - 1 to plot posterior draws and autocorrelations
%------------------------
prefs.map_smth = 1;     % - 1 to MAP smoothing (SLOW), 
prefs.ps_smth = 1;      % - 1 to perform particle smoothing (VERY SLOW),
                        %   both smoothers are O(nop^3) algorithms hence
                        %   are feasible only when "nop" is SMALL (<= 99)  
%%%%%%%%%%%%%%%%%%%%%%%%%

% data loading and preprocessing
[x,xm,x_opt,y,ym,y_opt,sta] = fn_proc_filen_(filename,prefs);

% processing x coordinate
[X_fin,X_diff,X_names] = fn_proc_cord_(filename,'X',x,xm,x_opt,del,nop,...
    resx,prefs);

% processing y coordinate
[Y_fin,Y_diff,Y_names] = fn_proc_cord_(filename,'Y',y,ym,y_opt,del,nop,...
    resy,prefs);

% save all results
eval(['save ',filename,'_results x xm x_opt X_fin X_diff ',...
    'X_names y ym y_opt Y_fin Y_diff Y_names sta filename del nop ',...
    'resx resy prefs'])

% this is a function which serves to embedd coordinates on frames used 
% and can be used after proper location of mat file with frames used
% is specified
%path2frames = 'D:\MPIBR\frames_used';
%fn_compare_(path2frames,filename,sta,x_opt,xm,X_fin,y_opt,ym,Y_fin);