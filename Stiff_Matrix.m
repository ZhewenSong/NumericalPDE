Err = [];
for n = 10:10:100; % # of elements
h = ones(n,1)/n; % element size, can be non-uniform mesh
grid = 10; % grid # in each element
x = zeros(n+1,1); 
x(1)=0; x(n+1)=1;
for i=2:n
    x(i)=x(i-1)+h(i);
end
A = zeros(2*(n-1),2*(n-1)); % initializing Stiff Matrix
A(1:2,1:2) = [12 -6*h(1);
           -6*h(1) 4*h(1)^2]/h(1)^3;    
for i=2:n-1
    j=2*i-3;
    a = [12 6*h(i) -12 6*h(i);
    6*h(i) 4*h(i)^2 -6*h(i) 2*h(i)^2;
    -12 -6*h(i) 12 -6*h(i);
    6*h(i) 2*h(i)^2 -6*h(i) 4*h(i)^2]/h(i)^3;
    A(j:j+3,j:j+3) = A(j:j+3,j:j+3) + a; 
end

A(2*n-3:2*n-2,2*n-3:2*n-2) = A(2*n-3:2*n-2,2*n-3:2*n-2) + ...
[ 12 6*h(n);
6*h(n) 4*h(n)^2]/h(n)^3;

F = zeros(2*(n-1),1); 
f = @(x) exp(x).*(x.^4 + 14*x.^3 + 49*x.^2 + 32*x - 12); % load function

ele_old = linspace(x(1),x(2),grid);

% initializing basis functions
phi0 = zeros(n-1,2*grid-1); 
phi1 = zeros(n-1,2*grid-1);

C_old = [1 x(1) x(1)^2 x(1)^3
        0 1 2*x(1) 3*x(1)^2
        1 x(2) x(2)^2 x(2)^3
        0 1 2*x(2) 3*x(2)^2];
CB0=C_old\[0;0;1;0];
CB1=C_old\[0;0;0;1];

phiB0_old = CB0(1)+CB0(2)*ele_old+CB0(3)*ele_old.^2+CB0(4)*ele_old.^3;
phiB1_old = CB1(1)+CB1(2)*ele_old+CB1(3)*ele_old.^2+CB1(4)*ele_old.^3;

% Calculating the loading vector based on the basis functions
for i=2:n
    C_new=[1 x(i) x(i)^2 x(i)^3
        0 1 2*x(i) 3*x(i)^2
        1 x(i+1) x(i+1)^2 x(i+1)^3
        0 1 2*x(i+1) 3*x(i+1)^2];
    
    % Coefficients of the cubic functions
    CA0=C_new\[1;0;0;0];
    CA1=C_new\[0;1;0;0];
    CB0=C_new\[0;0;1;0];
    CB1=C_new\[0;0;0;1];
    
    ele_new = linspace(x(i),x(i+1),grid);
    % Constructing the basis functions
    phiA0_new = CA0(1)+CA0(2)*ele_new+CA0(3)*ele_new.^2+CA0(4)*ele_new.^3;
    phiA1_new = CA1(1)+CA1(2)*ele_new+CA1(3)*ele_new.^2+CA1(4)*ele_new.^3;
    phiB0_new = CB0(1)+CB0(2)*ele_new+CB0(3)*ele_new.^2+CB0(4)*ele_new.^3;
    phiB1_new = CB1(1)+CB1(2)*ele_new+CB1(3)*ele_new.^2+CB1(4)*ele_new.^3;
    
    phi0(i-1,1:grid) = phiB0_old;
    phi0(i-1,grid:2*grid-1) = phiA0_new;
    phi1(i-1,1:grid) = phiB1_old;
    phi1(i-1,grid:2*grid-1) = phiA1_new;
    
    % Integration
    F(2*i-3)= sum(f([ele_old(1:grid-1),ele_new]).*phi0(i-1,:))*h(i)/(grid-1);
    F(2*i-2)= sum(f([ele_old(1:grid-1),ele_new]).*phi1(i-1,:))*h(i)/(grid-1);
    
    phiB0_old = phiB0_new;
    phiB1_old = phiB1_new;
    ele_old = ele_new;

end

u = A\F; % Calculating the coefficients of u
u0 = u(1:2:2*n-3);
u1 = u(2:2:2*n-2);
U = zeros(n,grid);
X = zeros(n,grid);

% Constructing solution u
for i=1:n-1
    part = u0(i)*phi0(i,:)+u1(i)*phi1(i,:);
    X(i,:) = linspace(x(i),x(i+1),grid);
    U(i,:) = U(i,:) + part(1:grid);
    U(i+1,:) = U(i+1,:) + part(grid:2*grid-1);
end
X(n,:) = linspace(x(n),x(n+1),grid);

plot(x,[U(:,1);0],'r-','linewidth',3); hold on
plot(x,exp(x).*x.^2.*(1-x).^2,'b:','linewidth',3)
legend('FEM soln','Exact soln')
%err = [U(:,1);0] - exp(x).*x.^2.*(1-x).^2;
err = zeros(n,grid);
for i=1:n
	err(i,:) = abs(U(i,:)-exp(X(i,:)).*X(i,:).^2.*(1-X(i,:)).^2);
end
Err = [Err log10(max(max(err)))];
end