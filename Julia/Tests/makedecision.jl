function makedecision(p_r,p_l,p_rope)
#resdec
#wl wr sl sr e
resdec=zeros(1,5)
if p_r>0.90
resdec[1,4]=1
elseif p_l>0.9
resdec[1,3]=1
elseif p_rope>0.95
resdec[1,5]=1
elseif p_r<0.05
resdec[1,1]=1
elseif p_l<0.05
resdec[1,2]=1

end

return resdec

end
