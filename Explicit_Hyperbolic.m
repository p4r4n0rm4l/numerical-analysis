% Õðïåñãáóßá 6
clear;
clc;

xmax = 1;
h = 0.1;
k = 0.1;
c = 1;
a = c*k/h;
t = 0.5;
x = 0:h:xmax;
n = round(xmax/h)+1;
nn = round(t/k)+1;
u = zeros(nn,n);

u(1,:) = 1/8*sin(pi*x); % Áñ÷éêÝò ÓõíèÞêåò
V = 0*x; 

for i=2:n-1
    u(2,i) = a^2*(u(1,i-1)+u(1,i+1))/2+(1-a^2)*u(1,i)+k*V(i);
end
for j=2:nn-1
    for i=2:n-1
        u(j+1,i) = a^2*(u(j,i-1)+u(j,i+1))+2*(1-a^2)*u(j,i)-u(j-1,i);
    end
end

surf(u);  % Ôï mesh ôçò u

D = 1/8*sin(pi*x)*cos(pi*t); % ÁíáëõôéêÞ ëýóç
E = abs(D'-u(nn,:)'); % Error
disp('   ÁíáëõôéêÞ:          Õðïëïãéóèåßóá:      ÓöÜëìá:   ');
disp([u(nn,:)' D' E]);
