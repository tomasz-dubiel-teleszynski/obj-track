function dsr = fn_ds_(parameters,y,sm,sv,kfr)

% Durbin-Koopman disturbance smoother

[T,Q,Z,H] = fn_sms(parameters,sm);

obs = size(y,1);
m = size(T,1);
p = size(y,2);

vt = kfr.yhat';
F1t = kfr.S1;
Kt = kfr.K;

ut = zeros(p,obs);
ust = zeros(p,obs);
epst = zeros(p,obs);

ksit = zeros(m,obs);

Dt{obs,1} = []; % p x p 
Vepst{obs,1} = []; % p x p 
Vksit{obs,1} = []; % m x m

rt = zeros(m,obs+1);
Nt{obs+1,1} = [];

rt(:,obs+1) = zeros(m,1);
Nt{obs+1} = zeros(m);

for t = obs:-1:sv.t0+1
    
    ut(:,t) = F1t{t} * vt(:,t) - Kt{t}'*T' * rt(:,t+1);
    epst(:,t) = H * ut(:,t);
    ksit(:,t) = Q * rt(:,t+1); 
    
    Dt{t} = F1t{t} + Kt{t}' *T'* Nt{t+1} * T* Kt{t};
    Vepst{t} = H - H * Dt{t} * H;
    Vksit{t} = Q - Q * Nt{t+1} *Q;   
    
    ust(:,t) = chol(inv(Dt{t}),'lower') * ut(:,t); 
    
    rt(:,t) = Z' * ut(:,t) + T' * rt(:,t+1);
    Nt{t} = Z' * Dt{t} * Z + T' * Nt{t+1} * T - Z' * Kt{t}'* T'* Nt{t+1} * T - T' * Nt{t+1} * T* Kt{t} * Z;

end

% output
dsr.ut = ut;
dsr.ust = ust;
dsr.epst = epst;
dsr.ksit = ksit;
dsr.Dt = Dt;
dsr.Vepst = Vepst;
dsr.Vksit = Vksit;
dsr.rt = rt;
dsr.Nt = Nt;

end
%% fn_sms
function [F,Q,H,R] = fn_sms(parameters,sm)

F = sm.A;
sig2 = exp(parameters(1));
Q = eval(sm.C(sig2));
H = sm.H;
expr = 'R = eval(sm.R(';
n = 1;
expr = strcat(expr,num2str(exp(parameters(n+1))));
for n = 2:size(H,1)
    expr = strcat(expr,strcat(',',num2str(exp(parameters(n+1)))));
end
expr = strcat(expr,'));');
eval(expr);

end