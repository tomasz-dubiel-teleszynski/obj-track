function pfr = fn_pf_(parameters,y,sm,sv,nop,res)

fprintf('\nParticle filtering begins...\n')

% create system matrices --------------------------------------------------
[F,Q,H,R] = fn_sms(parameters,sm);

T = size(y,1);
m = size(F,1);
S = zeros(T,m);
logL = zeros(T,1);

% used in smoothing
X1MH{T,1} = [];
P1MH{T,1} = [];
W1MH{T,1} = [];

% fill initial state values for coordinates -------------------------------
S = fn_fS(sv.t0,y,S,m);

% create particles --------------------------------------------------------
[X0,P0] = fn_cp(sv,H,m,nop,res); 

% create weights
[W0,w0] = fn_cw(nop);

for t = sv.t0+1:T
    
    %disp(t)
    
    if t < sv.t0+2
        X1 = X0;
        P1 = P0;
    end      
    
    % update particles ------------------------------------------------
    [X1,P1] = fn_up(y(t,:),X1,P1,F,Q,H,R,m);
    
    % calculate log likelihood for each particle ----------------------
    logLp = fn_logLp(y(t,:),X1,H,R,nop,res);
    
    % update and normalize weights
    W1 = fn_unw(w0,logLp);
    
    % calculate effective sample size
    ess = fn_ess(W1);
    
    if ess < nop/2;
        
        % resample particles ----------------------------------------------
        [X1,P1] = fn_rp(X1,P1,logLp,nop);
        
        if nop > 99 % rule of thumb
            
            % smoothing MH step -------------------------------------------
            [X1,P1] = fn_sMHs(y(t,:),X1,P1,F,Q,H,R,m,nop);
            
        end
        
        % reset weights
        W1 = W0;
        
    end
        
    % results -------------------------------------------------------------
    % state estimate
    S(t,:) = sum(kron(ones(m,1),W1).*X1,2);
    % likelihood estimate
    logL(t,1) = mean(logLp);
    % stored for smoothing
    X1MH{t} = X1;
    P1MH{t} = P1;
    W1MH{t} = W1;

end

% output
pfr.S = S;
pfr.logL = sum(logL);
pfr.X1MH = X1MH;
pfr.P1MH = P1MH;
pfr.W1MH = W1MH;

fprintf('\nParticle filtering done!\n')

end
%% fn_ess
function ess = fn_ess(W1)

ess = 1/sum(W1.^2);

end
%% fn_cw
function [W0,w0] = fn_cw(nop)

w0 = 1/nop*ones(1,nop);
W0 = w0;

end
%% fn_unw
function W1 = fn_unw(w0,logLp)

w1 = w0.*exp(logLp)';
W1 = w1/sum(w1);

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
%% fn_logLp
function logLp = fn_logLp(yt,X,H,R,nop,res)

logLp = zeros(nop,1);

for i = 1:nop
    Xi = X(:,i);
    C = H*Xi;
    if any(C>0 & C<=res)
        logLp(i,1) = fn_ll(yt,Xi,H,R);
    else
        logLp(i,1) = -Inf;
    end
end

end
%% fn_up
function [X,P] = fn_up(yt,X0,P0,F,Q,H,R,m)

X = zeros(size(X0));
P{size(X0,2),1} = [];
for i = 1:size(X0,2)
    [X(:,i),P{i}] = fn_upp(yt,X0(:,i),P0{i},F,Q,H,R,m);
end

end
%% fn_upp
function [X,P] = fn_upp(yt,X0,P0,F,Q,H,R,m)

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

X = mvnrnd(xtt,Ptt)';
P = Ptt;

end
%% fn_rp
function [X1,P1] = fn_rp(X,P,logLp,nop)

L = exp(logLp - max(logLp)); % likelihood
Q = L / sum(L); % normalization
R = cumsum(Q);
T = rand(1, nop);
[~, I] = histc(T, R);
idx=I+1;

X1 = X(:,I+1);
N = numel(P);
P1{N,1} = [];
for i = 1:N
    P1{i} = P{idx(i)};
end

end
%% fn_cp
function [X0,P0] = fn_cp(sv,H,m,nop,res)

X0 = zeros(m,nop);
P0{nop,1} = [];
for i = 1:nop
    X1 = mvnrnd(sv.x0,sv.P0)';
    C1 = H*X1;
    while ~all(C1>0 & C1<=res)
        X1 = mvnrnd(sv.x0,sv.P0)';
        C1 = H*X1;
    end
    X0(:,i) =  X1;
    P0{i} = sv.P0;
end

end
%% fn_ll
function ll = fn_ll(yt,Xp,H,R)

k = size(H,1);
C = H*Xp;
ll = - 0.5*log((2*pi)^k) - 0.5*log(det(R)) - 0.5*((yt'-C)'*inv(R)*(yt'-C));

end
%% fn_sMHs
function [X2,P2] = fn_sMHs(yt,X,P,F,Q,H,R,m,nop)

X2 = zeros(size(X));
P2{size(X,2),1} = [];
for i = 1:nop
    
    X0 = X(:,i);
    P0 = P{i}; 
    [X1,P1] = fn_upp(yt,X0,P0,F,Q,H,R,m);
    
    l0 = fn_ll(yt,X0,H,R);
    l1 = fn_ll(yt,X1,H,R);
    if log(rand()) <= l1-l0
        X2(:,i) = X1;
        P2{i} = P1; 
    else
        X2(:,i) = X0;
        P2{i} = P0;
    end
end 

end