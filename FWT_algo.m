function s = FWT_algo(v)
# Calculates the forward fast wavelet transform (FWT) of the given vector
# of single-scale coefficients, with respect to the biorthogonal b-spline 
# wavelet basis with d = d_tilde = 2

pkg load signal
c = v;
J = log2(length(c));
s = [];

a_filter = [-0.25, 0.5, 1.5, 0.5, -0.25];
b_filter = [0.5, -1, 0.5]; 

for i=1:J
  ci = (2^-0.5) * cconv(c, a_filter, 2^(J-i+1))(1:2:2^(J-i+1));
  ci = circshift(ci,-1);
  di = (2^-0.5) * cconv(c, b_filter, 2^(J-i+1))(1:2:2^(J-i+1));
  di = circshift(di,-1);
  c = ci;
  s = [di s];
endfor 

s = [ci s];
endfunction  