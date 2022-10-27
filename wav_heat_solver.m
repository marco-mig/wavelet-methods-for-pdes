function c_sol = wav_heat_solver(J)
# Calculates an approximation to the 1-D heat equation 
# u_t - mu * u_xx = f(x,t) 
# Using a multiscale wavelet basis (piecewise linear) corresponding to the 
# Biorthogonal B-spline Wavelets for d = d_tilda = 2  
# Based on the methods in the book "Scientific Computing with MATLAB and Octave"
# by Quarteroni (2014) 
# J corresponds to the level of the approximation 
# Example done for the given RHS f 

clear all
close all

# First we consider all the input parameters 
mu = 1; 
a = 0;
b = 1;
t1 = 0;
t2 = 1;
dt = 0.1;
f = @(x,t) -sin(2*pi*x)*sin(t) + (2*pi)^2*sin(2*pi*x)*cos(t);
theta = 1/2; 


# Value of u at t=0
u_t0 = @(x) sin(2*pi*x);

# Values of u at x=a and x=b (boundary points) 
u_0 = @(t) 0;
u_1 = @(t) 0; 


# We set up our wavelet elements 

fn{1} = @(x) phi(x); 
for j = 0:J
  for k = 0:2^j - 1
    fn{2^j+k+1} = @(x) 2^(j/2)*psi((2^j)*x-k);
  endfor
endfor 

for i=1:length(fn)
  fn{i}= @(x) fn{i}(x) + fn{i}(x-1) + fn{i}(x+1); 
endfor 

function y = phidv(x)
  ind1 = x >= -1 & x < 0;
  ind2 = x >= 0 & x < 1;
  ind3 = x < -1;
  ind4 = x >= 1;
  y(ind1) = 1; y(ind2) = -1; y(ind3) = 0; y(ind4) = 0;
endfunction 

function y = psidv(x)
  ind1 = x < -1;
  ind2 = x >= -1 & x < 0;
  ind3 = x >= 0 & x < 1/2;
  ind4 = x >= 1/2 & x < 1;
  ind5 = x >= 1 & x < 2; 
  ind6 = x >= 2;
  y(ind1) = 0; y(ind2) = -1/2; y(ind3) = 4; y(ind4) = -4; y(ind5) = 1/2; y(ind6) = 0;
endfunction 

dv{1} = @(x) phidv(x);
for j = 0:J
  for k = 0:2^j - 1 
    dv{2^j+k+1} = @(x) 2^(3*j/2)*psidv((2^j)*x-k);
  endfor
endfor 

for i=1:length(dv)
  dv{i} = @(x) dv{i}(x-2) + dv{i}(x-1) + dv{i}(x) + dv{i}(x+1) + dv{i}(x+2);
endfor 
  
N = length(fn) ;


# Calculating the mass matrix M using Gaussian Quadrature
M = zeros(N);
temp = @(x) (phi(x))^2;
M(1,1) = (1/2)*(temp((1/2)*(-3^(1/2)/3+1)) + temp((1/2)*(3^(1/2)/3+1)));
for j=0:J
  for k=0:2^j-1 
    for l=0:J
      for m=0:2^l-1 
        a1= (2^-j)*(k-1);
        b1= (2^-j)*(k+2);
        a2= (2^-l)*(m-1);
        b2= (2^-l)*(m+2);
        if b1-1 > 0
          c1 = 0;
        elseif b2-1 > 0
          c1 = 0;
        else
          c1 = max(0, min(a1,a2));
        endif
        
        if a1+1 < 1
          c2 = 1;
        elseif a2+1 < 1
          c2 = 1;
        else
          c2 = min(1, max(b1,b2));
        endif
        # First we find the matrix values when both functions are wavelets 
        
        temp = @(x) fn{2^j+k+1}(x) * fn{2^l+m+1}(x); 
        int = 0;
        q = 2^-(max(j,l)+1);
        grid = c1:q:c2; 
        for i=1:length(grid)-1 
          d1 = grid(i);
          d2 = grid(i)+q;
          int = int + (d2-d1)/2 * (temp((1/2)*((-(3^(1/2))/3)*(d2-d1) + d2 + d1)) + temp((1/2)*(((3^(1/2))/3)*(d2-d1) + d2 + d1)));
        endfor
        
        M(2^j+k+1, 2^l+m+1) = int;
      endfor
    endfor
  endfor
endfor

# We consider separately the cases with the scaling function       
for j=0:J
  for k=0:2^j-1
    temp = @(x) fn{2^j+k+1}(x) * fn{1}(x); 
    q = 2^(-j-1);
    grid = 0:q:1;
    int = 0;
    for i=1:length(grid)-1
      d1 = grid(i);
      d2 = grid(i)+q;
      int = int + (d2-d1)/2 * (temp((1/2)*((-(3^(1/2))/3)*(d2-d1) + d2 + d1)) + temp((1/2)*(((3^(1/2))/3)*(d2-d1) + d2 + d1)));
    endfor
    M(1, 2^j+k+1) = int;
    M(2^j+k+1, 1) = int;
  endfor
endfor  

# Calculating the "stiffness" matrix A 

A = zeros(N);
temp = @(x) (phidv(x))^2;
A(1,1) = mu*(1/2)*(temp((1/2)*(-3^(1/2)/3+1)) + temp((1/2)*(3^(1/2)/3+1)));
for j=0:J
  for k=0:2^j-1 
    for l=0:J
      for m=0:2^l-1 
        a1= (2^-j)*(k-1);
        b1= (2^-j)*(k+2);
        a2= (2^-l)*(m-1);
        b2= (2^-l)*(m+2);
        if b1-1 > 0
          c1 = 0;
        elseif b2-1 > 0
          c1 = 0;
        else
          c1 = max(0, min(a1,a2));
        endif
        
        if a1+1 < 1
          c2 = 1;
        elseif a2+1 < 1
          c2 = 1;
        else
          c2 = min(1, max(b1,b2));
        endif
        # First we find the matrix values when both functions are wavelets 
        
        temp = @(x) dv{2^j+k+1}(x) * dv{2^l+m+1}(x); 
        int = 0;
        q = 2^-(max(j,l)+1);
        grid = c1:q:c2; 
        for i=1:length(grid)-1 
          d1 = grid(i);
          d2 = grid(i)+q;
          int = int + (d2-d1)/2 * (temp((1/2)*((-(3^(1/2))/3)*(d2-d1) + d2 + d1)) + temp((1/2)*(((3^(1/2))/3)*(d2-d1) + d2 + d1)));
        endfor
        
        A(2^j+k+1, 2^l+m+1) = mu*int;
      endfor
    endfor
  endfor
endfor

# We consider separately the cases with the scaling function       
for j=0:J
  for k=0:2^j-1
    temp = @(x) fn{2^j+k+1}(x) * fn{1}(x); 
    q = 2^(-j-1);
    grid = 0:q:1;
    intb = 0;
    for i=1:length(grid)-1
      d1 = grid(i);
      d2 = grid(i)+q;
      int = int + (d2-d1)/2 * (temp((1/2)*((-(3^(1/2))/3)*(d2-d1) + d2 + d1)) + temp((1/2)*(((3^(1/2))/3)*(d2-d1) + d2 + d1)));
    endfor
    
    A(1, 2^j+k+1) = int;
    A(2^j+k+1, 1) = int;
  endfor
endfor  

# Now on to time discretization

t_grid = t1:dt:t2; 
h = 2^-(J+1);
x_grid = a:h:b;  

# Create the solution matrix u 
u = zeros(N, length(t_grid));

# We fill in the values that we know: 
# We have u_t0 (so we know t=0), aka the first column 

rhs = zeros(N,1);
for i=1:N
  rhs(i) = quad(@(x) u_t0(x)*fn{i}(x), a, b);
endfor 

col_1 = M\rhs;
u_k = col_1;
u(:,1) = u_k;

# We know u_0(t) = 0 and u_1(t) = 0, so the first and last rows are all zero 
# (i.e. Dirichlet boundary conditions) 
#u(1,:) = zeros(1,length(t_grid));
#u(length(x_grid),:) = zeros(1,length(t_grid));

f_n = zeros(N,1);
f_n1 = zeros(N,1);

for j=1:N
  f_n(j) = quad(@(x) f(x, t_grid(1))*fn{j}(x), a, b);
endfor

# Now we calculate the unknown values 

for i=2:length(t_grid)
  for j=1:N
  f_n1(j) = quad(@(x) f(x,t_grid(i))*fn{j}(x), a, b);  
  endfor 

  A1 = (1/dt)*M + theta*A;
  c1 = theta*f_n1 + (1-theta)*f_n - A*(1-theta)*u_k + (1/dt)*M*u_k;
  u(:,i) = A1 \ c1;
  
  u_k = u(:,i);
  f_n = f_n1;
endfor 

# Calculates u(x,t) at the initial time t=t1
initial = @(x) 0;
for i=1:N
  initial = @(x) initial(x) + u(i,1)*fn{i}(x); 
endfor 

# Calculates the approximation of u(x,t) at the final time t=t2
c_sol = @(x) 0;
for i=1:N
  c_sol = @(x) c_sol(x) + u(i, length(t_grid))*fn{i}(x);
endfor 

# The theoretical solution u(x,t) 
solution = @(x,t) sin(2*pi*x)*cos(t);

xi=0:2^-(J+1):1;
xi_fine = 0:0.0001:1;

subplot(1,2,1)
plot(xi_fine, solution(xi_fine,0),'--', xi, initial(xi), '*-')
title('Solution at t=0')
subplot(1,2,2)
plot(xi_fine, solution(xi_fine,1), '--', xi, c_sol(xi), '*-')
title('Solution at t=1')

endfunction 



