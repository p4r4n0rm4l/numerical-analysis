clear;
clc;

h = 0.2;
k = 0.01;
r = k/h^2;
xmax = 2;
tmax = 0.25;
n = floor(xmax/h)+1;
nn = floor(tmax/k)+1;

u = zeros(n,nn);
C = zeros(n-2,n-2);
B = zeros(n-2,n-2);

for i=1:n
    u(i,1) = sin(2*pi*(i-1)*h);
end

for j=1:nn
    u(1,j) = 0;
    u(n,1) = 0;
end

for i=1:n-2
    C(i,i) = 1+r;
    B(i,i) = 1-r;
end

for i=1:n-3
    C(i,i+1) = -r/2;
    C(i+1,i) = -r/2;
    B(i,i+1) = r/2;
    B(i+1,i) = r/2;
end

for j=1:nn-1
    w(2,j) = r*u(1,j);
    w(n-1,j) = r*u(n,j);
end

for j=1:nn-1
    u(2:n-1,j+1) = C\(B*(u(2:n-1,j)+w(2:n-1,j)));
end

surfl(u)
