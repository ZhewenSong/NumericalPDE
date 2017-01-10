% Defining A matrix of the linear equation system
A=[1 -1 1/2 -1/6 1/24
1 -1/2 1/8 -1/48 1/384
1 0 0 0 0
1 1 1/2 1/6 1/24
1 2 2 4/3 2/3]';
eh = [];
for h=[1/4,1/8,1/16,1/32]
    % solving for the coefficients c_j
    C = A\[0 0 1 0 0]'/h^2; 
    % error of numerical vs analytical
    eh = [eh, abs(C(1)*sin(2*(1-h))+C(2)*sin(2*(1-h/2))... 
    +C(3)*sin(2)+C(4)*sin(2*(1+h))+C(5)*sin(2*(1+2*h))+4*sin(2))]; 
end
scatter(log([1/4,1/8,1/16,1/32]),log(eh));
xlabel('log(h)')
ylabel('log(e_h)')
