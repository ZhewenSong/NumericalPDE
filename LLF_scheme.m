function LLF_scheme(U1,U2,x,J,k,h)
    % defining f(U)
    f1 = @(u1,u2) u2;
    f2 = @(u1,u2) u2.^2./u1+u1;

    for n=1:0.25/k
        for j=2:J-1
            %%% sup|f'(U)| for [U_j,U_{j+1}]
            ap=max([abs(U2(j:j+1)./U1(j:j+1)+1),abs(U2(j:j+1)./U1(j:j+1)-1)]);

            %%% Flux F_{j+1/2}
            Fp1=(f1(U1(j+1),U2(j+1))+f1(U1(j),U2(j)))/2-ap/2*(U1(j+1)-U1(j));
            Fp2=(f2(U1(j+1),U2(j+1))+f2(U1(j),U2(j)))/2-ap/2*(U2(j+1)-U2(j));

            %%% sup|f'(U)| for [U_{j-1},U_j]
            am=max([abs(U2(j-1:j)./U1(j-1:j)+1),abs(U2(j-1:j)./U1(j-1:j)-1)]);

            %%% Flux F_{j-1/2}
            Fm1=(f1(U1(j),U2(j))+f1(U1(j-1),U2(j-1)))/2-am/2*(U1(j)-U1(j-1));
            Fm2=(f2(U1(j),U2(j))+f2(U1(j-1),U2(j-1)))/2-am/2*(U2(j)-U2(j-1));

            %%% Conservative law
            U1(j)=U1(j)-k/h*(Fp1-Fm1);
            U2(j)=U2(j)-k/h*(Fp2-Fm2);

        end
    end
    plot(x,U1,'b:',x,U2./U1,'r:','linewidth',3); hold on