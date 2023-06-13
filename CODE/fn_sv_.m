function sv = fn_sv_(y,del)

%%%%%%%%%
t0 = 3;%%
%%%%%%%%%

y1 = y(t0,:);
y11 = diff(y(t0-1:t0,:),1)/del;
y12 = min(diff(y(1:t0,:),2)/(del^2));

vary = var(y);
vary1 = var(diff(y,1)/del);
vary2 = min(var(diff(y,2)/(del^2)));

[T,N] = size(y);
K = 2*min(T,N)+1;
x0 = zeros(1,K);
for k = 1:K-1
    if mod(k,2) == 0
        x0(k) = y11(1);
        y11(1)=[];
    else
        x0(k) = y1(1);
        y1(1) = [];
    end
end
x0(K)=y12;

P0 = eye(K);
for k = 1:K-1
    if mod(k,2) == 0
        P0(k,k) = vary1(1);
        vary1(1) = [];
    else
        P0(k,k) = vary(1);
        vary(1) = [];
    end
end
P0(K,K) = vary2;

sv.t0 = t0;
sv.x0 = x0;
sv.P0 = P0;

end