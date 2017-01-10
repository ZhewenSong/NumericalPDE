function Roe_scheme(U1,U2,x,J,k,h)
    vbar = @(rho_l,rho_r,v_l,v_r)... % The rho-averaged velocity
        (sqrt(rho_l)*v_l+sqrt(rho_r)*v_r)/(sqrt(rho_l)+sqrt(rho_r));

    for n=1:0.25/k
        for j=2:J-1
            %%% Roe flux F_{j+1/2}
            vbarp=vbar(U1(j),U1(j+1),U2(j)/U1(j),U2(j+1)/U1(j+1));
            [Fp1,Fp2]=Roe_flux(U1(j),U1(j+1),U2(j),U2(j+1),vbarp);
            
            %%% Roe flux F_{j-1/2}
            vbarm=vbar(U1(j-1),U1(j),U2(j-1)/U1(j-1),U2(j)/U1(j));
            [Fm1,Fm2]=Roe_flux(U1(j-1),U1(j),U2(j-1),U2(j),vbarm);
            
            %%% Conservative law
            U1(j)=U1(j)-k/h*(Fp1-Fm1);
            U2(j)=U2(j)-k/h*(Fp2-Fm2);
        end

    end
    plot(x,U1,'b--',x,U2./U1,'r--','linewidth',3); hold on