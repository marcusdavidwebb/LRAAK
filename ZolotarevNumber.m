function Z = ZolotarevNumber(a,b,c,d,n)
% assume -Inf <= a < b < c < d <= Inf
gamma = cross_ratio(a,b,c,d);
rho = exp(pi^2 / (pi * ellipke(1-1/gamma)/ellipke(1/gamma)));
Z = 4*rho.^(-n);
tau = 1;
fac = ((1+ rho.^(-4*tau*n)) ./ (1+rho.^(2*n*(1-2*tau)))).^4;
while fac < 1 - 1e-15
    Z = Z .* fac;
    tau = tau + 1;
    fac = ((1+ rho.^(-4*tau*n)) ./ (1+rho.^(2*n*(1-2*tau)))).^4;
end
Z(n==0) = 1;
end