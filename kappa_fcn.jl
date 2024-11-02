function kappa_fcn(x,y)
# if sqrt((x-1)^2+(y-0)^2)<0.25
#kappa=1E6
#elseif  sqrt((x-1)^2+(y-1)^2)<0.25
#   kappa=1E6
#else
#   kappa=0.0
#end

#if sqrt((x-1)^2+(y-0.5)^2)<0.25
#   kappa=1E6
#else
#   kappa=0.0
#end
   kappa=0.0
   xcenters=[1.25, 0.58, 1.51, 0.07, 0.73, 0.43, 1.83, 1.89, 0.02, 1.43, 1.59, 1.28, 1.03, 1.81, 0.61]
   ycenters=[0.33, 0.68, 0.66, 0.9, 0.4, 0.3, 0.69, 0.96, 0.56, 0.36, 0.17, 0.04, 0.08, 0.35, 0.43]
   R=0.05
   n_holes=size(xcenters,1)
   for i in 1:n_holes
      if sqrt((x-xcenters[i])^2+(y-ycenters[i])^2)< R
         kappa=1E7
      end
   end
   return kappa   
end