function V = eurocall(S0,sigma,X,r,T)
% EUROCALL(S0,sigma,X,r,T) calculates price at t = 0
% S0 -> num/array of initial stock prices

t = 0; % change this to calculate at any other time

N = @(d) 0.5*erf(d/sqrt(2)) + 0.5;
d1 = (log(S0/X) + (r+0.5*(sigma^2))*(T-t))/(sigma*sqrt(T-t));
d2 = (log(S0/X) + (r-0.5*(sigma^2))*(T-t))/(sigma*sqrt(T-t));

V = S0.*N(d1) - X*N(d2)*exp(-r*(T-t)); % for call

end

