function Ah = nabla2( Ah, x, y, n )

    stencil = [x-2:x+2; y-2:y+2];
    coeff = [1 -4 6 -4 1];
    for dim=1:2
        for i=1:5
            if stencil(dim,i)==-1
                stencil(dim,i) = n-1; %3;
            end
            if stencil(dim,i)==0
                stencil(dim,i) = n; %2;
            end
            if stencil(dim,i)==n+1
                stencil(dim,i) = 1; %n-1;
            end
            if stencil(dim,i)==n+2
                stencil(dim,i) = 2; %n-2;
            end
        end
    end
            

    for i=1:5
        Ah((x-1)*n+y,(stencil(1,i)-1)*n+stencil(2,3)) = ...
            Ah((x-1)*n+y,(stencil(1,i)-1)*n+stencil(2,3)) + coeff(i);
        Ah((x-1)*n+y,(stencil(1,3)-1)*n+stencil(2,i)) = ...
            Ah((x-1)*n+y,(stencil(1,3)-1)*n+stencil(2,i)) + coeff(i);
    end

    
    
end


