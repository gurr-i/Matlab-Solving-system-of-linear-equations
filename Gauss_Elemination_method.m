%%Gauss Elemination
clc
clear
format long g
fprintf('Gauss Elemination\n');
info = [24 -2860 7.26*10.^5;-2860 7.26*10.^5 -1.86472*10.^8;...
    7.26*10.^5 -1.86472*10.^8 5.24357*10.^10];
b = [1.057 * 10^(-4);-1.04162 * 10.^(-2);2.56799];
A = [info b];
for i=1:size(A,1)
    for j=i+1:size(A,1)
        key1 = A(j,i) ./ A(i,i);
        A(j,:) = A(j,:) - key1 .* A(i,:);
    end
end
%disp(A);
x = zeros(1,size(info,2));
for i = size(A,1): -1 : 1
    HG = sum(A(i,i+1:end-1).*x(i+1:end));
    x(i) =(A(i,end) - HG) ./A(i,i);
    %disp(x);
end
fprintf('Solution:  %d \n',x)
