function Bsigntest(y,x,rope)


s=1
nsamples=100000
diff=y-x

if rope>0
#Compute counts
nright=length(find(z->z> rope,diff[:]))
nleft=length(find(z->(z<-rope),diff[:]))
nrope=length(find(z->(z<=rope && z>=-rope),diff[:]))

data = rand(Dirichlet([nleft+s/3 nrope+s/3 nright+s/3 ][:]),nsamples)

else
#Compute counts
nright=length(find(z->z> 0,diff[:]))
nleft=length(find(z->(z<=0),diff[:]))

if nright<nleft

data = rand(Dirichlet([nleft nright+s][:]),nsamples)
elseif nright>nleft
data = rand(Dirichlet([nleft+s nright][:]),nsamples)

else

data = rand(Dirichlet([nleft+s/2 nright+s/2][:]),nsamples)


end

end
return data

end
