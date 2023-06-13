function [x,newests] = fn_outl_adj_dums_(ests,x,sm,sv,dm)

fprintf('\nMaximum likelihood estimation with dummies begins...\n')

options = fn_mle_opts_();
parsd = fn_init_parsd_(ests,x,dm);

[estsds, logL, exitflag] = fminsearch(@fn_logL_d,parsd,options,...
    x,sm,sv,dm);

fprintf('\nMaximum likelihood estimation with dummies done!\n')
fprintf(['Outcome = ', num2str(exitflag), '\n'])
fprintf(['Log likelihood = ', num2str(logL), '\n'])

x = (x'-fn_dm_C_(estsds,dm)*dm.D')'; 
newests = estsds(1:end-size(dm.D,2));

end
%% fn_init_parsd_
function parsd = fn_init_parsd_(ests,x,dm)

beta = std(x(:)); % outlier size
parsd = [ests; beta*ones(size(dm.D,2),1) ]; % dummies filled

end
%% fn_dm_C_
function C = fn_dm_C_(newests,dm)

expr1 = 'C = eval(dm.C(';
nparsd = numel(newests);
n0 = nparsd-size(dm.D,2);
expr1 = strcat(expr1,num2str(newests(n0+1)));
for n = n0+1:nparsd-1
    expr1 = strcat(expr1,strcat(',',num2str(newests(n+1))));
end
expr1 = strcat(expr1,'));');
eval(expr1);

end