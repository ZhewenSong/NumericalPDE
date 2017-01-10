function coeff = IFFT(func,j)
    % func is an array of n;
    % j is the j^th coeff;
    % T is period;
    n=size(func,2);
    if n>1
        if ~mod(j,2)
            coeff=IFFT(func(1:n/2),j/2)+IFFT(func(n/2+1:n),j/2);            
        else
            coeff=IFFT(func(1:n/2).*exp(sqrt(-1)*2*pi/n.*(1:n/2)),(j-1)/2)...
            -IFFT(func(n/2+1:n).*exp(sqrt(-1)*2*pi/n.*(1:n/2)),(j-1)/2);
        end
    else
        coeff=func(1);
    end
end