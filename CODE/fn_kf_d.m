function kfr = fn_kf_d(parametersd,y,sm,sv,dm)

% Kalman filter
% x = F*x1 + B + w, w ~ N(0,Q)
% z = H*x + v, v ~ N(0,R)

[F,Q,H,R,C,D] = fn_smsd(parametersd,sm,dm);

T = size(y,1);
p = size(y,2);
m = size(F,1);

yhat = zeros(T,p);
ust = yhat;
S{T,1} = []; % p x p
S1 = S;
K{T,1} = []; % m x p
xtt = zeros(T,m);
xtt1 = xtt;
Ptt{T,1} = []; % m x m
Ptt1 = Ptt;
logL = zeros(T,1);

% fill initial state values for coordinates -------------------------------
xtt = fn_fS(sv.t0,y,xtt,m);

% start values
xtt1(sv.t0+1,:) = sv.x0;
Ptt1{sv.t0+1,1} = sv.P0;

% filtering
for t = sv.t0+1:T
    
    if t > sv.t0+1
        
        % predicted state estimate
        xtt1(t,:) = xtt(t-1,:)*F';
        
        % predicted covariance
        Ptt1{t,1} = F*Ptt{t-1,1}*F' + Q;
    
    end
    
    % measurement residual
    yhat(t,:) = y(t,:) - xtt1(t,:)*H' - D(t,:)*C'; 
    
    % innovation covariance
    S{t,1} = H*Ptt1{t,1}*H' + R;
    S1{t,1} = inv(S{t,1});
    
    ust(t,:) = (sqrt(S1{t,1})*yhat(t,:)')';
    
    % Kalman gain
    K{t,1} = Ptt1{t,1}*H' *S1{t,1};
    
    % updated state estimate
    xtt(t,:) = xtt1(t,:) + yhat(t,:)*K{t,1}';
    
    % updated estimate covariance
    Ptt{t,1} = Ptt1{t,1} - K{t,1}*H*Ptt1{t,1};
    
    % log likelihood
    logL(t,1) = -0.5*( p*log(2*pi) + log(det(S{t,1})) + yhat(t,:)*S1{t,1}*yhat(t,:)' );
    
end
logL = sum(logL);

kfr.yhat = yhat;
kfr.S = S;
kfr.S1 = S1;
kfr.K = K;
kfr.xtt = xtt;
kfr.xtt1 = xtt1;
kfr.Ptt = Ptt;
kfr.Ptt1 = Ptt1;
kfr.logL = logL;
kfr.ust = ust;

end
%% fn_fS
function S = fn_fS(t0,y,S,m)

y1 = y(1:t0,:);

idx = logical([mod(1:m-1,2),0]);
for t = 1:t0
    S(t,idx) = y1(t,:);
end

end
%% fn_sms
function [F,Q,H,R,C,D] = fn_smsd(parametersd,sm,dm)

F = sm.A;
sig2 = exp(parametersd(1));
Q = eval(sm.C(sig2));
H = sm.H;
expr = 'R = eval(sm.R(';
n = 1;
expr = strcat(expr,num2str(exp(parametersd(n+1))));
for n = 2:size(H,1)
    expr = strcat(expr,strcat(',',num2str(exp(parametersd(n+1)))));
end
expr = strcat(expr,'));');
eval(expr);

D = dm.D;
nparametersd = numel(parametersd);
n0 = nparametersd-size(D,2);

expr1 = 'C = eval(dm.C(';
n = n0;
expr1 = strcat(expr1,num2str(parametersd(n+1)));
for n = n0+1:nparametersd-1
    expr1 = strcat(expr1,strcat(',',num2str(parametersd(n+1))));
end
expr1 = strcat(expr1,'));');
eval(expr1);

end