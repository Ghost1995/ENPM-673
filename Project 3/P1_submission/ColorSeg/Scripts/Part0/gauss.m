function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma)

    N = (1/(2*pi)^(length(x)/2))*(1/sqrt(det(sigma)))*exp(-0.5*((x - mu)/sigma)*(x - mu)');
    if isnan(N)
        N = 0;
    end
    
end