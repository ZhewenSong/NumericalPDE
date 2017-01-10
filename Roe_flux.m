function [F1,F2] = Roe_flux(u1_l,u1_r,u2_l,u2_r,vbar)
    
    % defining f(U)
    f1 = @(u1,u2) u2;
    f2 = @(u1,u2) u2.^2./u1+u1;
    
    % eigenvalues
    lambda1=vbar-1; 
    lambda2=vbar+1;
    
    % eigenvectors
    alpha1=1*(u1_r-u1_l)+(vbar-1)*(u2_r-u2_l);
    alpha2=1*(u1_r-u1_l)+(vbar+1)*(u2_r-u2_l);

    % Roe flux F=(f(u_l)+f(u_r))/2 - sum_p |lambda_p| alpha_p r_p
    F1=(f1(u1_l,u2_l)+f1(u1_r,u2_r))/2 ...
        -(abs(lambda1)*alpha1*1+abs(lambda2)*alpha2*1)/2;
    F2=(f2(u1_l,u2_l)+f2(u1_r,u2_r))/2 ...
        -(abs(lambda1)*alpha1*(vbar-1)+abs(lambda2)*alpha2*(vbar+1))/2;


end

