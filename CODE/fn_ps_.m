function psr = fn_ps_(pfr,y,parameters,sm,sv)

fprintf('\nParticle smoothing begins...\n')

[F,Q,H,R] = fn_sms(parameters,sm);

X1MH = pfr.X1MH;
P1MH = pfr.P1MH;
W1MH = pfr.W1MH;
Ss = pfr.S;

[T, m] = size(pfr.S);
nop = size(pfr.X1MH{sv.t0+1},2);

W{T,1} = [];

for t = sv.t0+1:T
    W{t} = W1MH{t};
end

WT{T,1} = W{T};

for t = T-1:-1:sv.t0+1
    %disp(t) 
    denom = zeros(nop,nop);
    for j = 1:nop
        for k = 1:nop
            denom(j,k) = W{t}(k)*exp( fn_lld(X1MH{t+1}(:,j),X1MH{t}(:,k),P1MH{t}{k},y(t,:),F,Q,H,R,m) );
        end
    end
    denom = sum(denom,2);
    for i = 1:nop
        brack =zeros(nop,1);
        for j = 1:nop
            brack(j) = WT{t+1}(j)/denom(j)*exp( fn_lld(X1MH{t+1}(:,j),X1MH{t}(:,i),P1MH{t}{i},y(t,:),F,Q,H,R,m) );
        end
        brack = sum(brack);
        WT{t}(i) = W{t}(i)*brack;
    end
    Ss(t,:) = sum(kron(ones(m,1),WT{t}).*X1MH{t},2);
end

psr.Ss = Ss;

fprintf('\nParticle smoothing done!\n')

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