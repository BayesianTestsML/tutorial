function heaviside(X)
#HEAVISIDE    Step function.
#    HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and .5 for X == 0.



Y = zeros(size(X));
Y[find( X-> X> 0,X)] = 1;
Y[find(X-> X== 0,X)] = .5;

return Y

end
