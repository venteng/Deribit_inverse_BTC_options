% filename: calc_inverse.m
% written by Huei-Wen Teng
% reviewed on 2022/8/31
% input
%   model
%   n = Monte Carlo sample size
%   U_base ~ Uniform([2*n, DTM_max]);
%   Z_base ~ Uniform([3*n, DTM_max]);


function price = calc_inverse(model, n, U_base, Z_base, omega, index_price, strike, r, y, d, maturity, param)

switch model
    
    case 'BS'
        
        price = calc_BS_inverse(omega, index_price, strike, r, y, d, maturity, param);
        
    case 'SV'
        
        mu = param(1);
        rho = param(2);
        alpha = param(3);
        beta = param(4);
        V0 = param(5);
        sig_v = param(6);
        
        Z = Z_base((1:(2*n)),:); % only need Z1, Z2
        
        S = get_SV(Z, param, index_price);
        payoff = calc_payoff_inverse(S, d, maturity, omega, strike,  r, index_price);
        price = mean(payoff); % model price
        
    case 'SVJ'
        
        
        % 1. define parameters from param
        mu = param(1);
        rho = param(2);
        alpha = param(3);
        beta = param(4);
        V0 = param(5);
        sig_v = param(6);
        lambda = param(7);
        mu_y = param(8);
        sig_y = param(9);
        
        U = U_base(1:n, :);   % Only need U, Z1, Z2, Z3
        Z = Z_base;
        
        S = get_SVJ(U, Z, param, index_price);
        payoff = calc_payoff_inverse(S, d, maturity, omega, strike,  r, index_price);
        price = mean(payoff); % model price
        
    case 'SVCJ'
        
        mu = param(1);
        rho = param(2);
        alpha = param(3);
        beta = param(4);
        V0 = param(5);
        sig_v = param(6);
        lambda = param(7);
        mu_y = param(8);
        rho_j = param(9);
        sig_y = param(10);
        mu_v = param(11);
        
        U = U_base;             % Need U1, U2
        Z = Z_base;             % Need Z1, Z2, Z3. 
        
        S = get_SVCJ(U, Z, param, index_price); % common random variable
        payoff = calc_payoff_inverse(S, d, maturity, omega, strike,  r, index_price);
        price = mean(payoff); % model price
        
    otherwise
        
        fprintf('warning in calc_inverse.m: No case!\n');
        price = 0;
        
end