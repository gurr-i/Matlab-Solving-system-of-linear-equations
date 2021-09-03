%%LU decomposition method for solving system of linear equation
clc
clear
close all
format long g
disp('LU decomposition');
% A = [1 2 5;0.2 1.6 7.4; 0.5 4 8.5];
A = [24 -2860 7.26.*10.^5;-2860 7.26.*10.^5 -1.86472.*10.^8;...
 7.26.*10.^5 -1.86472.*10.^8 5.24357.*10.^10];
B =[1.057*10^(-4);-1.04162*10^(-2);2.56799];
[L,U,P]=LU_decomposition_method(A);
y = L\(P*B);
x = U\y;
disp(x);

%%Function for LU_decomposition_method
function [L,U,P] = LU_decomposition_method(A)
NA = size(A,1);
AP = [A eye(NA)];
for k = 1:NA - 1
    [akx, kx] = max(abs(AP(k:NA,k)));
    if akx < eps
        error('Singular matrix and No LU decomposition')
    end
    mx = k+kx-1;
    if kx > 1 % Row change if necessary
        tmp_row = AP(k,:);
        AP(k,:) = AP(mx,:);
        AP(mx,:) = tmp_row;
    end
    % LU decomposition
    for m = k + 1: NA
        AP(m,k) = AP(m,k)/AP(k,k); % Eq.(2.4.8.2)
        AP(m,k+1:NA) = AP(m,k + 1:NA)-AP(m,k)*AP(k,k + 1:NA); % Eq.(2.4.9)
    end
end
P = AP(1:NA, NA + 1:NA + NA); % Permutation matrix
for m = 1:NA
    for n = 1:NA
        if m == n, L(m,m) = 1.; U(m,m) = AP(m,m);
        elseif m > n, L(m,n) = AP(m,n); U(m,n) = 0.;
        else L(m,n) = 0.; U(m,n) = AP(m,n);
        end
    end
end
if nargout == 0
    disp('L*U = P*A with');
    L,U,P
end
end
% You can check if P*L*U = A?
