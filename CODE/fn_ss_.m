function ssr = fn_ss_(kfr,y,sm,sv)

% Rauch–Tung–Striebel smoother
% xn = x + C*(x1n - x1) 
% Pn = P + C*(P1n - P1)*C'

F = sm.A;
T = size(y,1);

%m = size(F,1);

xtt = kfr.xtt;
xtt1 = kfr.xtt1;
Ptt = kfr.Ptt;
Ptt1 = kfr.Ptt1;

%----------
xtts = xtt;
Ptts = Ptt;
C = Ptts;
%----------

% % ---------------------------
% xtts(end,:) = xtt(end,:);
% Ptts{end,1} = Ptt{end,1};
% % ---------------------------

for t = T-1:-1:sv.t0+1
    
    C{t,1} = Ptt{t,1}*(F' /Ptt1{t+1,1});
    
    xtts(t,:) = xtt(t,:) + (xtts(t+1,:) - xtt1(t+1,:))*C{t,1}';
    Ptts{t,1} = Ptt{t,1} + C{t,1}*(Ptts{t+1,1} - Ptt1{t+1,1}*C{t,1})' ;
    
end

% output
ssr.xtts = xtts;
ssr.Ptts = Ptts;

end

