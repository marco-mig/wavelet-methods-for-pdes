function y = phi(x)
#PHI(X) Returns the biorthogonal B-Spline scaling function for (2,2) 
ind1= x<-1;
ind2= x>=-1 & x<0;
ind3= x>=0 & x<1;
ind4= x>=1; 
y(ind1)= 0;
y(ind2)= (x(ind2)+1);
y(ind3)= (-x(ind3)+1); 
y(ind4)= 0;
endfunction