clear;
clc;
format long;
an = @(x,y)(sin(5*pi*x/2)*sin(5*pi*y/2)); % Αναλυτική λύση
h = 0.05; 
nn = 0.4/h+1;   
num = zeros(nn);

k = 1;
for i=1:nn
    for j=1:nn
        num(i,j) = k;
        k = k+1;
    end
end

% Στον πίνακα n είναι οι κορυφές των τριγώνων 
n = zeros(2*(nn-1)^2,3);
k = 1;
for i=1:nn-1
    for j=1:nn-1
        n(k,1) = num(i,j);
        n(k,2) = num(i,j+1);
        n(k,3) = num(i+1,j);
        k = k+1;
        n(k,1) = num(i+1,j+1);
        n(k,2) = num(i+1,j);
        n(k,3) = num(i,j+1);
        k = k+1;
    end
end

% Ο coo έχει τις τιμές στου άξονες x,y
coo = zeros(nn^2,2);
k = 1;
for i=1:nn
    for j=1:nn
        coo(k,1) = (j-1)*h;
        coo(k,2) = (i-1)*h;
        k = k+1;
    end
end

% Ο πίνακας K έχει τις τιμές των M
KB = -(1/24)*(12.5*pi^2)*[2 1 1; 1 2 1; 1 1 2]*(h^2);
K = zeros(nn*nn);
no_of_elem = 2*4^(log2(nn-1));  % Αριθμός στοιχείων (τριγώνων)
for e=1:no_of_elem
    x(1) = coo(n(e,1),1);
    x(2) = coo(n(e,2),1);
    x(3) = coo(n(e,3),1);
    y(1) = coo(n(e,1),2);
    y(2) = coo(n(e,2),2);
    y(3) = coo(n(e,3),2);
    
    M(1,1) = ((y(2)-y(3))^2+(x(3)-x(2))^2)/(2*h^2);
    M(2,2) = ((y(3)-y(1))^2+(x(1)-x(3))^2)/(2*h^2);
    M(3,3) = ((y(1)-y(2))^2+(x(2)-x(1))^2)/(2*h^2);
    
    M(1,2) = (((y(2)-y(3))*(y(3)-y(1)))+(x(3)-x(2))*(x(1)-x(3)))/(2*h^2);
    M(2,1) = M(1,2);
    
    M(1,3) = (((y(2)-y(3))*(y(1)-y(2)))+(x(3)-x(2))*(x(2)-x(1)))/(2*h^2);
    M(3,1) = M(1,3);
    
    M(2,3) = (((y(3)-y(1))*(y(1)-y(2)))+(x(1)-x(3))*(x(2)-x(1)))/(2*h^2);
    M(3,2) = M(2,3);
    
    ke = -M+KB;
    for i=1:3
        for j=1:3
            K(n(e,i),n(e,j)) = K(n(e,i),n(e,j)) + ke(i,j);
        end
    end
end

%Βρίσκουμε τα boundaries και τα αφαιρούμε από τον K
bound = [num(1,:) num(2:(length(num)-1),1)' num(2:(length(num)-1),nn)' num(nn,:)];
K(bound,:) = [];
K(:,bound) = [];

% Ο πίνκας b έχει τις τιμές των p(i)
b = zeros(nn^2,1);
ksimax = @(k) (1-k);
const = -(h^2)*50*pi^2;
for e=1:no_of_elem
    x(1) = coo(n(e,1),1);
    x(2) = coo(n(e,2),1);
    x(3) = coo(n(e,3),1);
    y(1) = coo(n(e,1),2);
    y(2) = coo(n(e,2),2);
    y(3) = coo(n(e,3),2);

    g = @(l,k) (1-l-k).*(sin(2.5*pi*(x(1)+(x(2)-x(1)).*l+(x(3)-x(1)).*k)).*sin(2.5*pi*(y(1)+(y(2)-y(1)).*l+(y(3)-y(1)).*k)));
    be(1) = const*quad2d(g,0,1,0,ksimax);
    
    g = @(l,k) (l)    .*(sin(2.5*pi*(x(1)+(x(2)-x(1)).*l+(x(3)-x(1)).*k)).*sin(2.5*pi*(y(1)+(y(2)-y(1)).*l+(y(3)-y(1)).*k)));
    be(2) = const*quad2d(g,0,1,0,ksimax);
    
    g = @(l,k) (k)    .*(sin(2.5*pi*(x(1)+(x(2)-x(1)).*l+(x(3)-x(1)).*k)).*sin(2.5*pi*(y(1)+(y(2)-y(1)).*l+(y(3)-y(1)).*k)));
    be(3) = const*quad2d(g,0,1,0,ksimax);
    
    for i=1:3
        b(n(e,i)) = b(n(e,i))+be(i);
    end
    
end
 
% αφαιρούμε τα boundaries και από τον πίνακα b
b(bound) = [];

sol = K\b; % Λύση του συστήματος 

%mesh(reshape(sol,nn-2,nn-2)); % το mesh της λύσης

% Εύρεση αναλυτικής λύσης
ex = zeros(nn-2,nn-2);
for i=1:nn-2
    for j=1:nn-2
        ex(i,j)=an(i*h,j*h);
    end
end

err = abs(reshape(sol,nn-2,nn-2)-ex); % Εύρεση του σφάλματος υπολογισθείσας και αναλυτικής

%figure
%mesh(abs(reshape(sol,nn-2,nn-2)-ex)); % το mesh του σφάλματος