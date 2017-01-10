function Eel = stress(n, Rho)
    h = 1/(n-1);
    A = spalloc(2*n^2,2*n^2,18*n^2); 
    b = zeros(2*n^2,1);  
    e0 = 1;
    Rho = reshape(Rho,[n,n]);
    interp = tanh(Rho);
    %dinterp = sech(Rho).^2;
    for I=1:n
    for J=1:n
        A = deriv2(A,100*interp,1,1,1,1,I,J,n);
        A = deriv2(A,100*interp,2,2,2,2,I,J,n);
        A = deriv2(A,40*interp,1,1,2,2,I,J,n);
        A = deriv2(A,40*interp,2,2,1,1,I,J,n);
        A = deriv2(A,30*interp,1,2,1,2,I,J,n);
        A = deriv2(A,30*interp,2,1,2,1,I,J,n);
        A = deriv2(A,30*interp,1,2,2,1,I,J,n);
        A = deriv2(A,30*interp,2,1,1,2,I,J,n);
      
    end
    end


    for I=1:n
    for J=1:n
        Im = I-1; Ip = I+1;
        Jm = J-1; Jp = J+1;
        if I==1
            Im = n;
        elseif I==n
            Ip = 1;
        end
        if J==1
            Jm = n;
        elseif J==n
            Jp = 1;
        end
    b((I-1)*n+J) = (100+40)*(interp(Ip,J)-interp(Im,J))/2/h*e0*Rho(I,J);
    b(n^2+(I-1)*n+J)=(100+40)*(interp(I,Jp)-interp(I,Jm))/2/h*e0*Rho(I,J);
    end
    end

    U = A\b;
    U = reshape(U,[n,n,2]);
    
    Eel = zeros(n,n);
    for I=1:n
    for J=1:n
        Im = I-1; Ip = I+1;
        Jm = J-1; Jp = J+1;
        if I==1
            Im = n;
        elseif I==n
            Ip = 1;
        end
        if J==1
            Jm = n;
        elseif J==n
            Jp = 1;
        end
    Eel(I,J)=e0/2/h*(1.0*tanh(10*Rho(Ip,J)) + 0.4*tanh(10*Rho(Ip,J))) * ...
        (U(Ip,J,1)-U(Im,J,1) + U(I,Jp,2)-U(I,Jm,2));
    end
    end
    Eel = reshape(Eel,[n^2,1]);
end