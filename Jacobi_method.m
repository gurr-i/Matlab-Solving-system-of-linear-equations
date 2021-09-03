%%Jacobi method for solving system of linear equations
clc
clear
format long g
A = [24 -2860 7.26.*10.^5;-2860 7.26.*10.^5 -1.86472.*10.^8;...
 7.26.*10.^5 -1.86472.*10.^8 5.24357.*10.^10];
b =[1.057.*10.^(-4);-1.04162.*10.^(-2);2.56799];
fprintf('Jacobi\n');
n = length(b);
maxit = 1000;
x0 = rand(n,1);
tol = 1e-13;
x = Jacobi(A,b,tol,maxit,x0)
fprintf('MATLAB Built-in Function\n');
% making the augmented matrix
Awiggle = [A b];
temp = rref(Awiggle);
x = temp(:,end)

  
%%Function for Jacobi
function x = Jacobi(A, b, epsilon, maxit, x)
% check if the entered matrix is a square matrix
[na, ma] = size(A);
if na ~= ma
 disp('ERROR: Matrix A must be a square matrix')
 return
end
% check if b is a column matrix
[nb, mb] = size (b);
if nb ~= na || mb~=1
 disp( 'ERROR: Matrix b must be a column matrix')
 return
end
dx = zeros(na,1);
for k=1:maxit
 sum = 0;
 for i=1:na
  dx(i) = b(i);
 for j=1:na
  dx(i) = dx(i) - A(i,j)*x(j);
 end
 dx(i) = dx(i)/A(i,i);
 x(i) = x(i) + dx(i);
 if (dx(i) >= 0)
   sum = sum + dx(i);
 else
  sum = sum - dx(i);
 end
 end
 if(sum <= epsilon)
   break
 end
end
fprintf('The final answer obtained after %g iterations is \n', k);
end
