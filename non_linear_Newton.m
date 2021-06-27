% ������� Newton ��� �� �������� ���������
clear;
clc;
format long;

% �����������
F = inline('[x(1)+cos(x(1)*x(2)*x(3))-1; (1-x(1))^0.25+x(2)+0.05*x(3)^2-0.15*x(3)-1; -x(1)^2-0.1*x(2)^2+0.01*x(2)+x(3)-1]');

% ������� ���������
DF = inline('[1-x(2)*x(3)*sin(x(1)*x(2)*x(3)),-x(1)*x(3)*sin(x(1)*x(2)*x(3)),-x(1)*x(2)*sin(x(1)*x(2)*x(3)); -0.25*(1-x(1))^-0.75 , 1 , -0.1*x(3)-0.15; -2*x(1),-0.2*x(2)+0.01,1]');

% ������
tol = 1e-5;

% �������� ������� �����������
max_iter = 50;

% ������ ����������
x0 = [0.1;0.1;0.1];

i = 0;

while (i<=max_iter)
    i = i + 1;
    increment = -DF(x0)\F(x0);
    x = x0 + increment;
    diff = norm(x-x0,inf);
    disp([i x' diff]);
    if (diff<tol)
        break;
    end
    x0 = x;
end