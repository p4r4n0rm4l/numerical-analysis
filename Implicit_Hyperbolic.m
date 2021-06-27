clear;
clc;

L = 1;
h = 0.1;
k = 0.1;
c = 1;
a = c*k/h;
t = 0.5;
x = 0:h:L;

n = round(L/h)+1;
nn = round(t/k)+1;
u = zeros(n,nn);
u(:,1) = 1/8*sin(pi*x);
V = zeros(n,1);
V(1,:) = 0;
A = gallery('tridiag',n-2);
I = speye(n-2);
b = zeros(n,1);

b(1) = -1;
b(n) = -1;
u(2:n-1,2) = (2*(4/(a^2)*I+A))\((2*(4/(a^2)*I-A))*u(2:n-1,1)+(4/a^2*I+A)*2*k*V(2:n-1,1));
for j=2:nn-1
    u(2:n-1,j+1) = ((4/(a^2)*I+A))\((2*(4/(a^2)*I-A))*u(2:n-1,j)-(4/a^2*I+A)*u(2:n-1,j-1));
end
surf(u');

D = 1/8*sin(pi*x)*cos(pi*t);
E = abs(D'-u(:,nn))
