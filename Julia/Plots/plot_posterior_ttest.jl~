function plot_posterior_ttest(i,j,d,mur,sigmar,dofr,marginleft,marginright)

xs1=marginleft:0.0001:marginright
xs2=-0.01:0.001:0.01

f=x -> pdf(TDist(dofr[1]), (x-mur[1])/sigmar[1])

df1 = DataFrame(DeltaAcc=xs1,pdf=f(xs1),ymin=xs1*0,ymax=f(xs1),legend="pdf")
df2 = DataFrame(DeltaAcc=xs2,pdf=0.0001,ymin=xs2*0,ymax=0.0001,legend="rope")
df = vcat(df1, df2)

p=plot(df, x=:DeltaAcc, y=:pdf, ymin=:ymin, ymax=:ymax, color=:legend, xintercept=[-0.01, 0.01],Geom.vline(color="orange", size=1mm),Geom.line, Geom.ribbon,Theme(major_label_font_size=13pt,minor_label_font_size=12pt,key_label_font_size=11pt),Coord.Cartesian(xmin=marginleft,xmax=marginright))


draw(PDF("Plots/output$i$j$d.pdf", 6inch, 3inch), p)

return p

end
