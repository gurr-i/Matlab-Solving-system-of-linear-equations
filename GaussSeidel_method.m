%%GaussSeidel method for solving system of linear equations
clc
clear
format long g
disp('Gauss-Seidel');
A = [24 -2860 7.26.*10.^5;-2860 7.26.*10.^5 -1.86472.*10.^8;...
    7.26.*10.^5 -1.86472.*10.^8 5.24357.*10.^10];
b =[1.057.*10.^(-4);-1.04162.*10.^(-2);2.56799];

n = length(b);
x0 = rand(n,1);
tol = 1e-13;
x = Gauss_Siedel(A,b,x0,tol)
  
  
%%Function for Gauss Siedel Method
function x = Gauss_Siedel(A,b,x0,tol)
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
% Separation of matrix A into lower triangular and upper triangular matrices
% A = D + L + U
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
% check for convergence condition for Gauss-Seidel method
e= max(eig(-inv(D+L)*(U)));
if abs(e) >= 1
 disp (['Since the modulus of the largest Eigen value of ' ...
 'iterative matrix is not less than 1'])
 disp ('this process is not convergent.')
 return
end
k = 1;
x(:,1) = x0;
err = 1000000000*rand(na,1);% initial error assumption for looping
while sum(abs(err) >= tol) ~= zeros(na,1)
 x(:,k+1) = -inv(D+L)*(U)*x(:,k) + inv(D+L)*b;% Gauss-Seidel formula
 err = x(:,k+1) - x(:, k);% finding error
 k = k + 1;
end
fprintf('The final answer obtained after %g iterations is \n', k);
x = x(:,end);
end
