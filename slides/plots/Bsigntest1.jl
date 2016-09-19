function Bsigntest(y,x,rope,priorvec=[1 1 1]/3)


s=1
nsamples=100000
diff=y-x

if rope>0
#Compute counts
nright=length(find(z->z> rope,diff[:]))
nleft=length(find(z->(z<-rope),diff[:]))
nrope=length(find(z->(z<=rope && z>=-rope),diff[:]))

data = rand(Dirichlet([nleft nrope nright][:]+priorvec*s),nsamples)

else
#Compute counts
nright=length(find(z->z> 0,diff[:]))
nleft=length(find(z->(z<=0),diff[:]))


data = rand(Dirichlet([nleft nright][:]+priorvec*s),nsamples)



end
return data

end


