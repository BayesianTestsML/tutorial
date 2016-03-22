function Bttest_correlated(x,rho,m,rope_min_value,rope_max_value,prob)
#Bayesian correlated t-test
#x: differences of accuracies
#rho: correlation
#m: mean to be tested (m=0)
#the rope interval is [rope_min_value,rope_max_value]
#prob: probability for HDI credible interval

#Return parameters of the posterior tildemu, sigma_student, dof
# probability of the three regions right, left, rope p_right, p_left, p_rope,
# hdi interval

#parameters of the non-informative prior, to match the results of the
#frequentist correlated t-test
k0=1000000;
a=-0.5;
b=0;
mu0=0;
m=0;

#used for HDI computation
alpha=1-prob

if rope_min_value>rope_max_value
    error("rope_min_value larger than rope_max_value")
end


n=length(x)
H=ones(n,1)

#correlation matrix
M=rho*ones(n,n)+diagm(ones(n))*(1-rho)
#posterior mean
tildemu=inv((H'*inv(M)*H+1/k0))*(H'*inv(M)*x+mu0/k0)

#sample mean
hatmu=inv((H'*inv(M)*H))*(H'*inv(M)*x)

#posterior parameters of the Gamma distribution
tildea=a+n/2
tildeb0=0.5*((x-H*hatmu)'*inv(M)*(x-H*hatmu)+2*b-(H*hatmu)'*inv(M)*(H*hatmu)-mu0^2/k0+((H*tildemu)'*inv(M)*(H*tildemu)+tildemu'*tildemu/k0));
tildeb=max(tildeb0,0);

#parameters of the Student distribution
sigma_student=sqrt((tildeb/tildea/(H'*inv(M)*H+1/k0))); 
dof=2*tildea

#for numerical problems
if sigma_student[1]<0.000000001
	sigma_student[1]=0.000000001
end
t=(tildemu-m)/sigma_student;

#compute HDI interval for the posterior
crit =  invlogcdf(TDist(2*tildea),log((1 - alpha / 2))).* sigma_student;
hdi = [tildemu - crit; tildemu + crit];



#Computations of probability of three regions for accuracy=class1-class2: left (classif2 better than classif1)
# right (classif1 better than classif2) and rope (classif1=classif2)

p_left=0.5
p_right=0.5
p_rope=0.0;


p_right=1-cdf(TDist(2*tildea), (rope_max_value-tildemu)/sigma_student)
p_left=cdf(TDist(2*tildea), (rope_min_value-tildemu)/sigma_student)
p_rope=cdf(TDist(2*tildea), (rope_max_value-tildemu)/sigma_student)-cdf(TDist(2*tildea), (rope_min_value-tildemu)/sigma_student)



if rope_max_value==rope_min_value
    p_rope=0.0;
    p_left=1-p_right;
end



#return parameters of the posterior tildemu, sigma_student, dof
#probability of the three regions right, left, rope p_right, p_left, p_rope,
# hdi interval

return tildemu, sigma_student, dof, p_right, p_left, p_rope, hdi


end
