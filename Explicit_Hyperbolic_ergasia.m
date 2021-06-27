% Υποεργασία 6
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

u(1,:) = 1/8*sin(pi*x); % Αρχικές Συνθήκες
V = 0*x; 

for i=2:n-1
    u(2,i) = a^2*(u(1,i-1)+u(1,i+1))/2+(1-a^2)*u(1,i)+k*V(i);
end
for j=2:nn-1
    for i=2:n-1
        u(j+1,i) = a^2*(u(j,i-1)+u(j,i+1))+2*(1-a^2)*u(j,i)-u(j-1,i);
    end
end

surf(u);  % Το mesh της u

D = 1/8*sin(pi*x)*cos(pi*t); % Αναλυτική λύση
E = abs(D'-u(nn,:)'); % Error
disp('   Αναλυτική:          Υπολογισθείσα:      Σφάλμα:   ');
disp([u(nn,:)' D' E]);