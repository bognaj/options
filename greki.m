function [C, delta, gamma, theta, vega, rho] = greki(St, K, r, sigma, T, t, cp)

d_plus = (log(St/K) + (r + sigma^2/2)*(T-t))/(sigma*sqrt(T-t));
d_minus = (log(St/K) + (r - sigma^2/2)*(T-t))/(sigma*sqrt(T-t));

if cp == "call"
    C = St * normcdf(d_plus) - K * exp(-r*(T - t)) * normcdf(d_minus);
    delta = normcdf(d_plus);
    gamma = normpdf(d_plus)/(St*sigma*sqrt(T-t));
    theta = 1/(2*sqrt(T-t)) * St * normpdf(d_plus) * sigma + r * K * exp(-r*(T-t)) * normcdf(d_minus);
    vega = St * normpdf(d_plus) * sqrt(T - t);
    rho = K * (T - t) * exp(-r * (T-t)) * normcdf(d_minus);
    
elseif cp == "put"
    C = -St*normcdf(-d_plus) + K*exp(-r*(T-t)) * normcdf(-d_minus);
    delta = -normcdf(-d_plus);
    gamma = normpdf(d_plus)/(St*sigma*sqrt(T-t));
    theta = 1/(2*sqrt(T-t)) * St * normpdf(d_plus) * sigma - r * K * exp(-r*(T-t)) * normcdf(-d_minus);
    vega = St * normpdf(d_plus) * sqrt(T - t);
    rho = -K * (T - t) * exp(-r * (T-t)) * normcdf(-d_minus);
end
end