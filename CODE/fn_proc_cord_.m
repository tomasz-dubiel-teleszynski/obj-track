function [X_fin,X_diff,names] = fn_proc_cord_(filename,labl,x,xm,x_opt,del,nop,res,prefs)

fprintf('\nProcessing coordinate begins...\n');

% create system matrices;
sm = fn_sm_(x,del);

% initialize state - based on observations
sv = fn_sv_(x,del);

% maximum likelihood estimation
ests = fn_mle_(x,sm,sv);

% outlier detection
[dsr,dm] = fn_outl_det_(ests,x,sm,sv);

% outlier adjustment WITH MEANS (FAST)
xms = fn_outl_adj_mean_(x,dsr);

if isempty(prefs.adj_dums)
    prefs.adj_dums = input('\nIf outlier adjustment WITH DUMMIES (SLOW) is\nNOT TO BE RUN press 0, otherwise 1, and then ENTER!\n');
end
if prefs.adj_dums

    % outlier adjustment WITH DUMMIES (SLOW)
    xds = fn_outl_adj_dums_(ests,x,sm,sv,dm);

    % ourlier adjustment CHOICE (USER)
    newx = fn_outl_adj_choice_(x,xds,xms,prefs);
    
    % reinitialize state (1st) - based on observations
    newsv = fn_sv_(newx,del);
    
    % maximum likelihood reestimation - based on states
    newests = fn_mle_re_(newx,sm,newsv);
        
    % reinitialize state (2nd)
    kfr = fn_kf_(newests,newx,sm,newsv);
    idx = fn_idx_(newx);
    newsv = fn_sv_(kfr.xtt(:,logical(idx)),del);
    
else
    
    % automatic choice of outlier adjustment WITH MEANS (FAST)
    newx = xms;
    
    % reinitialize state (1nd) - based on observations
    newsv = fn_sv_(newx,del);
    
    % maximum likelihood reestimation
    newests = fn_mle_re_(newx,sm,newsv);
    
    % reinitialize state (2nd) - based on states
    kfr = fn_kf_(newests,newx,sm,newsv);
    idx = fn_idx_(newx);
    newsv = fn_sv_(kfr.xtt(:,logical(idx)),del);
    
end

if isempty(prefs.mcmc_ests)
    prefs.mcmc_ests = input('\nIf Bayesian MCMC is TO BE RUN press 1,\notherwise press 0, and then ENTER!\n');
end
if prefs.mcmc_ests
    
    % MCMC estimation
    finests = fn_mcmce_(newests,newx,del,prefs);
    
else
    
    % automatic choice of maximum likelihood estimates
    finests = newests;
    
end

% Kalman FILTERING (FAST)
kfr = fn_kf_(finests,newx,sm,newsv);

% Kalman SMOOTHING (FAST)
ssr = fn_ss_(kfr,newx,sm,newsv);

% particle FILTERING (FAST)
pfr = fn_pf_(finests,newx,sm,newsv,nop,res);

if isempty(prefs.map_smth)
    prefs.map_smth = input('\nIf MAP (SLOW) is NOT TO BE RUN press 0,\notherwise 1, and then ENTER!\n');
end
if prefs.map_smth

    % MAP (SLOW)
    map = fn_map_(pfr,newx,finests,sm,newsv);

end

if isempty(prefs.ps_smth)
    prefs.ps_smth = input('\nIf PARTICLE SMOOTHING (VERY SLOW) is NOT TO BE RUN press 0,\notherwise 1, and then ENTER!\n');
end
if prefs.ps_smth

    % particle SMOOTHING (VERY SLOW)
    psr = fn_ps_(pfr,newx,finests,sm,newsv);

end

% final coordinate estimates;
if prefs.ps_smth && prefs.map_smth
    
    names = [{'KF'},{'SS'},{'PF'},{'PS'},{'MAP'}];
    [X_fin,X_diff] = fn_fin_cord_(filename,labl,x,xm,x_opt,finests,names,kfr.xtt,ssr.xtts,pfr.S,psr.Ss,map.Ss);

elseif prefs.ps_smth && ~prefs.map_smth

    names = [{'KF'},{'SS'},{'PF'},{'PS'}];
    [X_fin,X_diff] = fn_fin_cord_(filename,labl,x,xm,x_opt,finests,names,kfr.xtt,ssr.xtts,pfr.S,psr.Ss);
    
elseif ~prefs.ps_smth && prefs.map_smth
    
    names = [{'KF'},{'SS'},{'PF'},{'MAP'}];
    [X_fin,X_diff] = fn_fin_cord_(filename,labl,x,xm,x_opt,finests,names,kfr.xtt,ssr.xtts,pfr.S,map.Ss);
    
else
    
    names = [{'KF'},{'SS'},{'PF'}];
    [X_fin,X_diff] = fn_fin_cord_(filename,labl,x,xm,x_opt,finests,names,kfr.xtt,ssr.xtts,pfr.S);
    
end

fprintf('\nProcessing coordinate finished!\n');
fprintf('\n-------------------------------\n')

end
%% fn_outl_adj_choice
function x = fn_outl_adj_choice_(x,xds,xms,prefs)

if isempty(prefs.adj_choice)
    
    figure
    plot(x)
    title('\bfCoordinates without adjustment')
    legd = legend('X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8);
    
    figure
    plot(xds)
    title('\bfCoordinates adjusted with dummies')
    legd = legend('X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8);
    
    figure
    plot(xms)
    title('\bfCoordinates adjusted with means')
    legd = legend('X-correll','Cascade','Foreground');
    set(legd,'FontAngle','italic','FontSize',8);
    
    prefs.adj_choice = input('\nWhich adjustment you choose? With dummies press 1,\nwith means press 0, then press ENTER!\n');
    
end

if prefs.adj_choice
    
    x = xds;
    
else
    
    x = xms;
    
end

end
%% fn_idx_
function idx = fn_idx_(x)

idx = zeros(2*size(x,2)+1,1);
for n = 1:numel(idx)-1
    if mod(n,2) ~= 0
        idx(n)=1;
    end
end

end