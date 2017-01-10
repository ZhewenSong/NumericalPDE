clear; clc;

figure
h = 1/99;
k = h/10;
epsilon = h/2;
x = 0:h:1;
y = 0:h:1;
n = size(x,2);
[X,Y]=meshgrid(x,y);
%load('random_Rho_100.mat');
Rho = sin(X*50*pi).*sin(Y*50*pi);

Rho = reshape(Rho,n^2,1);

axis equal
cbar = colorbar;
set(cbar,'ylim',[-1,1])
D = spdiags(Rho,0,n^2,n^2);
Ah1 = spalloc(n^2,n^2,5*n^2); 
for i=1:n
    for j=1:n
        Ah1 = nabla1(Ah1,i,j,n);
    end
end
Ah1 = sparse(Ah1)/h^2;
Ah2 = spalloc(n^2,n^2,7*n^2); 
for i=1:n
    for j=1:n
        Ah2 = nabla2(Ah2,i,j,n);
    end
end
Ah2 = sparse(Ah2)/h^4;
T = 100;
M = zeros(T,1);
E = zeros(T,1);
ax = gca;
ax.NextPlot = 'replaceChildren';
F(T) = struct('cdata',[],'colormap',[]);
for t=1:T

    left = speye(n^2) + k*epsilon^2*Ah1^2 - k*Ah1*D^2;
    right = Rho - k*Ah1*Rho;% - k*Ah1*stress(n,Rho);
    Rho = left\right;
    D = spdiags(Rho,0,n^2,n^2);
    rho = reshape(Rho,n,n);
    M(t) = sum(Rho)*h^2;

    pcolor(X,Y,rho); shading interp; 
    drawnow;

    for i=1:n
        for j=1:n
        if i==1 
            im = n; 
        else im = i; 
        end
        if i==n
            ip = 1;
        else
            ip = i;
        end
        if j==1
            jm = n;
        else
            jm = j;
        end
        if j==n
            jp = 1;
        else
            jp = j;
        end
        E(t) = E(t) + epsilon^2/2 * ...
        (((rho(ip,j)-rho(im,j))/2/h)^2+((rho(i,jp)-rho(i,jm))/2/h)^2) + ...
        rho(i,j)^4/4-rho(i,j)^2/2;
        end
    end
    E(t) = E(t)*h^2;
    
end
