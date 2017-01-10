function coeff = FFT(func,k)
    % func is an array of n;
    % k is the k^th coeff;
    % T is period;
    n=size(func,2);
    if n>1
        if ~mod(k,2)
            coeff=FFT(func(1:n/2),k/2)+FFT(func(n/2+1:n),k/2);            
        else
            coeff=FFT(func(1:n/2).*exp(-sqrt(-1)*2*pi/n.*(1:n/2)),(k-1)/2)...
            -FFT(func(n/2+1:n).*exp(-sqrt(-1)*2*pi/n.*(1:n/2)),(k-1)/2);
        end
    else
        coeff=func(1);
    end
end

