function [dsr,dm] = fn_outl_det_(parameters,x,sm,sv)

kfr = fn_kf_(parameters,x,sm,sv);
dsr = fn_ds_(parameters,x,sm,sv,kfr);
dm = fn_dm_(dsr);

end
%% fn_dm_
function dm = fn_dm_(dsr)

ust = dsr.ust';
[N,p] = size(ust);
idx= logical(abs(ust)>3.291); % it is t-test's critical value for 99.9%
                              % whereas 2.576 is for 99%
D = [];
dimd = size(D,2);
nd = zeros(1,p);

x = zeros(N,1);
for j = 1:p
    [loc, pts] = find(idx(:,j)==1);
    npts = numel(pts);
    if npts < 2
        d = x;
        d(loc) = 1;
        D = [D,d];
    else
        start = 0;
        for l = 1:npts
            if loc(l)>start
                k = 1;
                while idx(loc(l)+k,j)==1
                    k = k + 1;
                end
                if k < 2
                    d = x;
                    d(loc(l)) = 1;
                    D = [D,d];
                else
                    d = x;
                    d(loc(l):loc(l)+k-1)=1;
                    D = [D,d];
                    start = start + loc(l) + k;
                end
            end
        end
    end
    nd(1,j) = size(D,2) - dimd;
    dimd = size(D,2);
end

Nd = sum(nd);
expr = 'syms';
BETA{Nd,1} = [];
for n = 1:Nd
    BETA{n,1} = strcat(' beta',num2str(n));
    expr = strcat(expr,BETA{n,1});
end
eval(expr);

c(p,Nd) = eval(BETA{Nd,1});
for i = 1:p
    if i < 2
        for j = 1:nd(i)
            c(i,j) = eval(BETA{j,1});
        end
    else        
        for j = sum(nd(1:i-1))+1:sum(nd(1:i))
            c(i,j) = eval(BETA{j,1});
        end
    end
end

expr1 = 'C = symfun(c,[';
for n = 1:Nd
    expr1 = strcat(expr1,BETA{n,1});
end
expr1 = strcat(expr1,']);'); 
eval(expr1);

dm.C = C;
dm.D = D;

end

