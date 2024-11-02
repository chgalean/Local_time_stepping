function kappa_fcn(x,y)
 if sqrt((x-1)^2+(y-0)^2)<0.25
    kappa=1E6
 elseif  sqrt((x-1)^2+(y-1)^2)<0.25
    kappa=1E6
 else
    kappa=0.0
 end
 return kappa   
end