function options = fn_mle_opts_()

tol = 1e-5;
iter = 10000;

% 'iter' to show each iteration of the minimization, 'off' otherwise
displ = 'off';

options = optimset('Display',displ,'MaxFunEval',iter,'MaxIter',iter,...
    'TolX',tol,'TolFun',tol);

end

