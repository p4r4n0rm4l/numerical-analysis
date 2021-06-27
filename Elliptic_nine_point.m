% ���������� 4 ������� �
clear;
clc;
format long;

box = 0.4;      % ������� ���� � ��� �
h = 0.05;        % ���� ��� � ��� � �����
n = box/h - 1;  % ������� ���������� �������

f = @(x,y)(-25*pi^2*sin(2.5*pi*x)*sin(2.5*pi*y));   % � f(x,y) ���������
an=@(x,y)(sin(2.5*pi*x)*sin(2.5*pi*y));             % ��������� ���� ��� ��������


% ��������� ������ �
A = eye(n^2)*(-20-12.5*pi^2*6*h^2+2*h^2*12.5*pi^2); % ������� ������ ���������

num = -(h^2/2)*12.5*pi^2;    % ��������� ��� �� 5-point molecule tou �u
for i=1:n^2-1                % ������� ��� ��� ���� ���������
    A(i,i+1) = 4 + num;
    A(i+1,i) = 4 + num;
    if (mod(i,n)==0)
        A(i+1,i) = 0;
        A(i,i+1) = 0;
    end
end

for i=1:n^2-n                
    A(i,n+i) = 4 + num;
    A(n+i,i) = 4 + num;
end

for i=1:(n^2-n-1)
    A(i,n+1+i) = 1;
    A(n+i+1,i) = 1;
    if (mod(i,n)==0)
        A(i,n+i+1) = 0;
        A(n+i+1,i) = 0;
    end
end

for i=1:n^2-n+1
    A(i,n+i-1) = 1;
    A(n+i-1,i) = 1;
    if (mod(i-1,n)==0)
        A(i,n+i-1) = 0;
        A(n+i-1,i) = 0;
    end
end

k = 1;
s = zeros(n^2,1);   % ���������� ������ S
for i=1:n           % ������� ��� S ������ �� ��� ��������� �����
    for j=1:n
        s(k,1) = (h^2/2)*(8*f(i*h,j*h)+f((i-1)*h,j*h)+f((i+1)*h,j*h)+f(i*h,(j-1)*h)+f(i*h,(j+1)*h));    
        k = k+1;
    end
end

u = A\s;            % ���� ��� ����������

uu = zeros(n^2,1);
k = 1;
for i=1:n
    for j=1:n
        uu(k) = an(i*h,j*h);    % ������ ���������� ����� ��� ��� ��������
        k = k+1;
    end
end


disp('   ���������:          �������������:      ������:   ');
disp([u uu abs(u-uu)])

k = 1;
for i=1:n
    for j=1:n
        sah(i,j) = u(k);
        k=k+1;
    end
end