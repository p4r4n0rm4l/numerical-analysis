% ���������� 5 ������� �
clear;
clc;
format long;

h = 0.2;
k = 0.01;
r = k/h^2;

xmax = 2; % ������� ���� ��� x
tmax = 0.25; % ������� ���� ��� t

n = floor(xmax/h)+1;
nn = floor(tmax/k)+1;

D = zeros(1,n-2);  % �� ����� ��� ������ ���������
UD = zeros(1,n-3); % �� ����� ��� ��� ���������
LD = zeros(1,n-3); % �� ����� ��� ���� ���������
u = zeros(n,nn);  % ������� ��� �������� ��� ����

for i=1:n
    u(i,1) = sin(2*pi*(i-1)*h);
end

for j=1:nn
    u(1,j) = 0;
    u(n,j) = 0;
end

% ������� ������ ���������
for i=1:n-2
    D(1,i) = 1+2*r;
end

% ������� ��� ��� ���� ���������
for i=1:n-3
    UD(1,i) = -r;
    LD(1,i) = -r;
end

for j=1:nn-1
    w(2,j) = r*u(1,j);
    w(n-1,j) = r*u(n,j);
end

% � ���� ��� ����������
u = sparse_choleski_parabolic_1D(D,UD,LD,u,nn,n);
surfl(u)