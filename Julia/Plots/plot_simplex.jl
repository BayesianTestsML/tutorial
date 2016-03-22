function projectSimplex(points)
    
    #Project probabilities on the 3-simplex to a 2D triangle
    
    #N points are given as N x 3 array
    ss=size(points,1)
    # Convert points one at a time
    tripts = zeros(size(points,1),2)
    for idx=1:size(points,1)
        # Init to triangle centroid
        x = 1.0 / 2
        y = 1.0 / (2 * sqrt(3))
        # Vector 1 - bisect out of lower left vertex 
        p1 = points[idx, 3]
        x = x - (1.0 / sqrt(3)) * p1 * cos(pi / 6)
        y = y - (1.0 / sqrt(3)) * p1 * sin(pi / 6)
        # Vector 2 - bisect out of lower right vertex  
        p2 = points[idx, 1]  
        x = x + (1.0 / sqrt(3)) * p2 * cos(pi / 6)
        y = y - (1.0 / sqrt(3)) * p2 * sin(pi / 6)        
        # Vector 3 - bisect out of top vertex
        p3 = points[idx, 2]
        y = y + (1.0 / sqrt(3) * p3)
        tripts[idx,:] = [x,y]
    end
    return tripts
        
    end
    

function plot_simplex(data,name1,name2)
dataproj=projectSimplex(data')


vert=projectSimplex([1 0 0;0 1 0;0 0 1])

vert0=projectSimplex([0.3333 0.3333 0.3333; 0.5 0.5 0;0.5 0 0.5;0 0.5 0.5])


df = DataFrame(Pleft=dataproj[:,1][:], Prope=dataproj[:,2][:])


p=plot(df, x=:Pleft, y=:Prope, Guide.xlabel(""), Guide.ylabel(""),  Guide.annotation(
         compose(context(), text(-0.04, -0.09, name1 ))),Guide.annotation(compose(context(), text(0.88, -0.09, name2 ))),Guide.annotation(compose(context(), text(0.44, 0.9, "rope" ))), Geom.hexbin,Theme(
key_label_font_size=11pt),  Coord.Cartesian(xmin=-0.1,xmax=1.05,ymin=-0.1,ymax=1.05),Guide.xticks(ticks=nothing), Guide.yticks(ticks=nothing), Guide.annotation(
compose(context(), polygon([(vert[1,1], vert[1,2]), (vert[2,1], vert[2,2]), (vert[3,1], vert[3,2])]), fill(nothing),
       stroke(colorant"orange"))), Guide.annotation(
compose(context(), polygon([(vert0[1,1], vert0[1,2]),(vert0[2,1], vert0[2,2]),(vert0[1,1], vert0[1,2]),(vert0[3,1], vert0[3,2]),(vert0[1,1], vert0[1,2]),(vert0[4,1], vert0[4,2])]), fill(nothing),
       stroke(colorant"orange"),linewidth(0.125mm))))
#
return p

end
