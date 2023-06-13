function logL = fn_logL_(parameters,y,sm,sv)

kfr = fn_kf_(parameters,y,sm,sv);
logL = -kfr.logL;

end

