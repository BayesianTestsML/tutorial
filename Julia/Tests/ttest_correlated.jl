function ttest_correlated(x,m,rho,tail,alpha)
#Frequentist correlated t-test
#x: differences of accuracies
#m: mean to be tested (m=0)
#rho: correlation
#tail: tail=0 (two-tailed), tail=1 (right one-tailed test), tail=2 (left one-tailed test)
#alpha: alpha level

#return p-value and confidence interval

samplesize=length(x)

#define the correction
te_tr_ratio=rho/(1-rho)

#degrees of freedom
df = max(samplesize - 1,0);

#statistics
xmean = mean(x);
sdpop = std(x);
#for numerical stability
if sdpop ==0
    sdpop=10.0^-10;
end

#correlation correction of the standard error
ser = sdpop .* sqrt(1/samplesize+te_tr_ratio);

#t statistics
tval = (xmean - m) ./ ser;


# Compute the correct p-value for the test, and confidence intervals
if tail == 0 # two-tailed test
    p = 2 * cdf(TDist(df),-abs(tval))
    crit =  invlogcdf(TDist(df),log((1 - alpha / 2))).* ser;
    ci = [xmean - crit; xmean + crit];
elseif tail == 1 # right one-tailed test
    p = cdf(TDist(df),-tval)
        crit = invlogcdf(TDist(df),log((1 - alpha ))).* ser;  
        ci = [xmean - crit; Inf];
else tail == -1 # left one-tailed test
    p = cdf(TDist(df),tval)
         crit = invlogcdf(TDist(df),log((1 - alpha ))).* ser;   
  ci = [-Inf; xmean + crit];
end

#return p-value and confidence interval
return p, ci


end
