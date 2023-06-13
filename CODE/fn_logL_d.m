function logL = fn_logL_d(parametersd,y,sm,sv,dm)

kfr = fn_kf_d(parametersd,y,sm,sv,dm);
logL = -kfr.logL;

end

