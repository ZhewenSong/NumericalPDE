function coeff = RFT(func,k,T)
    coeff=0;
    n=size(func,2);
    for j=0:n-1
        coeff=coeff+func(k+1)*exp(-sqrt(-1)*k*j*T/n);
    end
    %coeff=coeff/n;

end

