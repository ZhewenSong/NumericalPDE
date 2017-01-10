
h = 0.1;
k = 1e-7;
x = 0:h:1;
n = size(x,2);
b = 1;
gamma = 0.001;
C = 0;
e = 1;
u = rand(size(x));
T = 2000;
M = zeros(T,1);
for t=1:T
    for j=1:n
        stencil = j-2:j+2;
        for index = 1:5
            if stencil(index)==-1
                stencil(index) = 3;
            end
            if stencil(index)==0
                stencil(index) = 2;
            end
            if stencil(index)==n+1
                stencil(index) = n-1;
            end
            if stencil(index)==n+2
                stencil(index) = n-2;
            end
        end
        u(j) = u(j) + k*(u(stencil(2))^3-2*u(stencil(3))^3+u(stencil(4))^3)/h^2 ...
            - b^2*(u(stencil(2))-2*u(stencil(3))+u(stencil(4)))/h^2 ...
            - gamma*(u(stencil(1))-4*u(stencil(2))+6*u(stencil(3))-4*u(stencil(4))+u(stencil(5)))/h^4;
            %- C*e*(u(J5)-2*u(J4)-2*u(J2)+u(J1))/2/h^3); 
    end
    M(t) = sum(u)*h;
    plot(1:t,M(1:t));
    %ylim([0,1]);
    pause(0.25); 
end