% filename: calc_BS_delta.m
function delta = calc_BS_delta(omega, index_price, strike, r, y, d, maturity, param)
% European call option formula under the Black-Scholes model is correct.
% https://goodcalculators.com/black-scholes-calculator/
% This is a good webpage to check if the Black-Scholes formula is correct
% Parameters are annualized
% re-parametrization
sig = param;
sig2 = sig*sig; % annualized
strikehat = 1/strike;
tau = (maturity-d); % on a daily basis. All parameters are on a daily basis. 

% 
d1 = log(index_price/strike)/(sig*sqrt(tau))+(r-y+0.5*sig)*sqrt(tau);
d2 = d1-sig*sqrt(tau);
d3 = d2-sig *sqrt(tau);

% Corrected by Tom
delta = omega * (exp((sig2)*tau) * (1/index_price^2)*strike*normcdf(omega*d3));

% The old formula from Alexandra
%delta = omega * (exp((sig2)*tau) * (1/index_price^2)*strike*normcdf(omega*d2));

