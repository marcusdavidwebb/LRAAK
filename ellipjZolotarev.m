function [roots,poles,sigma] = ellipjZolotarev(a,b,c,d,n)
    % assume -Inf <= a < b < c < d <= Inf
    M = cross_ratio(a,b,c,d);
    gam = -1 + 2*M + 2*sqrt(M)*sqrt(M - 1);

    % T maps (-1,-1/gam,1/gam,1) -> (a,b,c,d)
    if isinf(a)
        z = [-1/gam;1/gam;1];
        w = [b;c;d];
    else
        z = [-1;-1/gam;1/gam];
        w = [a;b;c];
    end 
    t = [z(1)*w(1)*(w(2) - w(3)) + z(2)*w(2)*(w(3) - w(1)) + z(3)*w(3)*(w(1) - w(2));
         z(1)*w(1)*(z(2)*w(3)-z(3)*w(2))+z(2)*w(2)*(z(3)*w(1)-z(1)*w(3))+z(3)*w(3)*(z(1)*w(2)-z(2)*w(1));
         w(1)*(z(3)-z(2)) + w(2)*(z(1)-z(3)) + w(3)*(z(2)-z(1));
         z(1)*w(1)*(z(2) - z(3)) + z(2)*w(2)*(z(3) - z(1)) + z(3)*w(3)*(z(1) - z(2))];
    T = @(z) (t(1)*z + t(2)) ./ (t(3)*z + t(4));

    % Following https://github.com/ajt60gaibb/freeLYAP/blob/master/iterative_solvers/getshifts_adi.m, 
    % we use asymptotic approxmation when gam is large.
    if gam < 1e7
        K = ellipke(1-(1/gam)^2);
        u = (1/2:n-1/2)*K/n;
        [~, ~, dn] = ellipj(u, 1-(1/gam)^2);
    else
        K = (2*log(2)+log(gam)) + (-1+2*log(2)+log(gam))/gam^2/4;
        u = (1/2:n-1/2)*K/n;
        dn = sech(u) + (0.25/gam^2)*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u);
    end

    roots = T(dn.');
    poles = T(-dn.');
    sigma = ZolotarevNumber(a,b,c,d,n);
end