function [estimates,logL,exitflag] = fn_mle_re_(x,sm,sv)

fprintf('\nMaximum likelihood reestimation begins...\n')

options = fn_mle_opts_();
pars = fn_init_pars(x);

[estimates, logL, exitflag] = fminsearch(@fn_logL_,pars,options,...
    x,sm,sv);

fprintf('\nMaximum likelihood reestimation done!\n')
fprintf(['Outcome = ', num2str(exitflag),'\n'])
fprintf(['Log likelihood = ', num2str(logL),'\n'])

end
%% fn_init_pars_
function pars = fn_init_pars(x)

pars(1,1) = 1; % sig2
pars(2:2+size(x,2)-1,1) = 1; % noise

end