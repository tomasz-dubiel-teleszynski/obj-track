function x = fn_outl_adj_mean_(x,dsr)

ust = dsr.ust';
idx= ~logical(abs(ust)>3.291); 

[n, m] = size(idx);

for i = 2:n-1
    for j = 1:m
        if idx(i,j) == 0
            k=1;
            while idx(i+k,j) == 0
                k = k+1;
            end
            x(i,j)= ( x(i-1,j) + x(i+k,j) ) / 2;
        end
    end
end

