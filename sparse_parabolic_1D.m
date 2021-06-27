% Υποεργασία 5 ερώτημα α
clear;
clc;
format long;

h = 0.2;
k = 0.01;
r = k/h^2;

xmax = 2; % Μέγιστη τιμή του x
tmax = 0.25; % Μέγιστη τιμή του t

n = floor(xmax/h)+1;
nn = floor(tmax/k)+1;

D = zeros(1,n-2);  % Οι τιμές της κύριας διαγωνίου
UD = zeros(1,n-3); % Οι τιμές της άνω διαγωνίου
LD = zeros(1,n-3); % Οι τιμές της κάτω διαγωνίου
u = zeros(n,nn);  % Πίνακας που περιέχει την λύση

for i=1:n
    u(i,1) = sin(2*pi*(i-1)*h);
end

for j=1:nn
    u(1,j) = 0;
    u(n,j) = 0;
end

% Γέμισμα κύριας διαγωνίου
for i=1:n-2
    D(1,i) = 1+2*r;
end

% Γέμισμα ’νω και Κάτω διαγωνίου
for i=1:n-3
    UD(1,i) = -r;
    LD(1,i) = -r;
end

for j=1:nn-1
    w(2,j) = r*u(1,j);
    w(n-1,j) = r*u(n,j);
end

% Η λύση του συστήματος
u = sparse_choleski_parabolic_1D(D,UD,LD,u,nn,n);
surfl(u)