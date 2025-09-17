function gamma = cross_ratio(a,b,c,d)
    % assume -Inf <= a < b < c < d <= Inf
    if isinf(a)
        gamma = abs((d-b)/(c-b));
    elseif isinf(d)
        gamma = abs((c-a)/(c-b));
    else
        gamma = abs(((c-a)*(d-b))/((c-b)*(d-a)));
    end 
end