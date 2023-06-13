function sm = fn_sm_(y,del)

[T,N] = size(y);
K = 2*min(T,N)+1;

syms arg_sig2 arg_noise2 t s d;

a(K,K) = 0;
for n = 1:K-1
    if mod(n,2) == 0
        a(n,end) = 1;
    else
        a(n,n+1) = 1;
    end
end
c(K,K) = 1;

A = expm(del*a);

% no need for "(c*c')" here below since "c" is enough 
tempC = int(expm((t-s)*a)*c*expm((t-s)*a'),'s','t-d','t');
tempC = symfun(tempC,d); 
C = symfun(arg_sig2*tempC(del),arg_sig2);

% OLD VERSION
% tempC = simplify(int(expm((t-s)*a)*c*expm((t-s)*a'),'s',strcat('t-',num2str(del)),'t'));
% C = symfun(arg_sig2*tempC,arg_sig2);

sm.A = A;
sm.C = C;

H = zeros(N,K); 
for n = 1:N
    temp = zeros(1,K);
    temp(2*n-1) = 1;
    H(n,:) = temp;
end

sm.H = H;

expr = 'syms';
NOISE{N,1} = [];
for n = 1:N
    NOISE{n,1} = strcat(' noise',num2str(n));
    expr = strcat(expr,NOISE{n,1});
end
eval(expr);
r(N,N) = arg_sig2; % arg_sig used here just to define symbolic matrix
for n = 1:N
    r(n,n) = NOISE{n,1};
end
expr1 = 'R = symfun(r,[';
for n = 1:N
    expr1 = strcat(expr1,NOISE{n,1});
end
expr1 = strcat(expr1,']);'); 
eval(expr1);

sm.R = R;

end