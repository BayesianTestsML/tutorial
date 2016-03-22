function get_accuracies(i,j,d,ClassID,DatasetID,Percent_correct)

indi=find(x->x==i,ClassID)
indj=find(x->x==j,ClassID)
indd=find(x->x==d,DatasetID)
indid=intersect(indi,indd)
indjd=intersect(indj,indd)
acci=Percent_correct[indid]/100
accj=Percent_correct[indjd]/100

return acci, accj

end
