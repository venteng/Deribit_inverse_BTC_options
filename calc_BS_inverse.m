% filename: calc_BS_inverse.m
% written by Huei-Wen Teng
% reviewed on 2022/8/31
function price = calc_BS_inverse(omega, index_price, strike, r, y, d, maturity, param)
sig = param;
sig2 = sig*sig; % annualized
strikehat = 1./strike; % could be a vector or scalor
tau = (maturity-d); % in daily. % all parameters are on a daily basis.

if tau > 0
    
    d1 = log(index_price./strike)/(sig*sqrt(tau))+(r-y+0.5*sig)*sqrt(tau); % could be a vector or scalor
    d2 = d1-sig * sqrt(tau); % could be a vector
    d3 = d2-sig * sqrt(tau); % could be a vector
    price = omega .* ( exp(-r*tau)*normcdf(omega.*d2) - exp((y-r+sig2)*tau) .* (1./index_price).*strike.*normcdf(omega.*d3));
else
    
    Shat = 1./index_price;
    price =  exp(-r*tau) * Shat .* max( omega .*(index_price - strike), 0);    % could be a vector or scalor
end