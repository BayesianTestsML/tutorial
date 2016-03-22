function plot_post_hierch(x,marginleft,marginright,title)




df = DataFrame(mu0=x)

#Geom.histogram(bincount=20
p=plot(df, x=:mu0,  xintercept=[-0.01, 0.01],Guide.title(title),Geom.vline(color="orange", size=1mm), Geom.density,Coord.Cartesian(xmin=marginleft,xmax=marginright),Theme(major_label_font_size=13pt,minor_label_font_size=12pt,key_label_font_size=11pt))



return p

end
