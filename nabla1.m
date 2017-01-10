function Ah = nabla1( Ah, x, y, n )

    stencil = [x-1:x+1; y-1:y+1];
    coeff = [1 -2 1];
    
    for dim=1:2
        for i=1:3
            
            if stencil(dim,i)==0
                stencil(dim,i) = n; 
            end
            if stencil(dim,i)==n+1
                stencil(dim,i) = 1; 
            end
            
        end
    end

    for i=1:3
        Ah((x-1)*n+y,(stencil(1,i)-1)*n+stencil(2,2)) = ...
            Ah((x-1)*n+y,(stencil(1,i)-1)*n+stencil(2,2)) + coeff(i);
        Ah((x-1)*n+y,(stencil(1,2)-1)*n+stencil(2,i)) = ...
            Ah((x-1)*n+y,(stencil(1,2)-1)*n+stencil(2,i)) + coeff(i);
    end

    
    
end