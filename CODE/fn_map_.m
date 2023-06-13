function psr = fn_map_(pfr,y,parameters,sm,sv)

fprintf('\nMAP smoothing begins...\n')

[F,Q,H,R] = fn_sms(parameters,sm);

X1MH = pfr.X1MH;
P1MH = pfr.P1MH;
Ss = pfr.S;

[T, m] = size(pfr.S);
nop = size(pfr.X1MH{sv.t0+1},2);

delta = zeros(T,nop);

t = sv.t0+1;
for i = 1:nop
    delta(t,i) = fn_llp(X1MH{t}(:,i),sv) + fn_ll(y(t,:),X1MH{t}(:,i),H,R);
end

theta = zeros(T,nop);

for t = sv.t0+2:T
    %disp(t)
    for j = 1:nop
        logp = zeros(1,nop);
        for i = 1:nop
             logp(i) = fn_lld(X1MH{t}(:,j),X1MH{t-1}(:,i),P1MH{t-1}{i},y(t,:),F,Q,H,R,m);
        end
        brack = delta(t-1,:) + logp;
        maxbrack = max(brack);
        delta(t,j) =  fn_ll(y(t,:),X1MH{t}(:,j),H,R) + maxbrack;
        theta(t,j) = find(brack==maxbrack,1,'first');
    end
end

idx = zeros(T,1);

% for t = T;
idx(T) = find(delta(T,:)==max(delta(T,:)),1,'first');
Ss(T,:) = X1MH{T}(:,idx(T));

for t = T-1:-1:sv.t0+1
    idx(t) = theta(t+1,idx(t+1));
    Ss(t,:) = X1MH{t}(:,idx(t));
end

psr.Ss = Ss;

fprintf('\nMAP smoothing done!\n')

end
%% fn_ll
function ll = fn_ll(yt,Xp,H,R)

k = size(H,1);
C = H*Xp;
ll = - 0.5*log((2*pi)^k) - 0.5*log(det(R)) - 0.5*((yt'-C)'*inv(R)*(yt'-C));

end
%% fn_llp
function llp = fn_llp(Xp,sv)

k = numel(Xp);
llp = - 0.5*log((2*pi)^k) - 0.5*log(det(sv.P0)) - 0.5*((Xp-sv.x0')'*inv(sv.P0)*(Xp-sv.x0'));

end
%% fn_lld
function lld = fn_lld(X1,X0,P0,yt,F,Q,H,R,m)

[xtt,Ptt] = fn_upp(yt,X0,P0,F,Q,H,R);

lld = - 0.5*log((2*pi)^m) - 0.5*log(det(Ptt)) - 0.5*((X1-xtt)'*inv(Ptt)*(X1-xtt));

end
%% fn_upp
function [xtt,Ptt] = fn_upp(yt,X0,P0,F,Q,H,R)

xtt1 = F*X0;
Ptt1 = F*P0*F' + Q;

S = H*Ptt1*H' + R;
K = Ptt1*(H'/S);
yhat = yt' - H*xtt1;    

xtt = xtt1 + K*yhat;
Ptt = Ptt1 - K*H*Ptt1;

% "safety net"--------------
Ptt = 0.5*(Ptt+Ptt');
if det(Ptt)==0
    Ptt = Ptt + 1e-6*eye(m);
end
% --------------------------

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