function y = psi(x) 
#PSI(X) returns the Biorthogonal B-Spline Wavelet (2,2)
ind1 = x<-1;
ind2 = x>=-1 & x<0;
ind3 = x>=0 & x<1/2;
ind4 = x>=1/2 & x<=1;
ind5 = x>1 & x<=2;
ind6 = x>2;
y(ind1)= 0;
y(ind2)= -(-x(ind2)/2 - 1/2); 
y(ind3)= -(4*x(ind3) - 1/2);
y(ind4)= -(-4*x(ind4) + 7/2);
y(ind5)= -(x(ind5)/2 - 1);
y(ind6)= 0;
endfunction