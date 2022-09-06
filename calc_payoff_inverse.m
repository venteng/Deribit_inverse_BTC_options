% filename: calc_payoff.m
% written by Huei-Wen Teng (2021/8/26, 2021/10/20)
% input
% cp: types of the Bitcoin option: {'C', 'P', 'EC', 'EP'}
% strike: a scalar, strike price
% d: a scalar, current day
% maturity: a scalar, maturity day
% r: a scalar, interest rate
% index_price: initial stock price to pervent the case of DTM = 0.
% S: a matrix historical stock prices path of size (n x DTM)
% output
% payoff: a vector of size (n x 1)

function payoff = calc_payoff_inverse(S, d, maturity, omega, strike,  r, index_price)

% reparametrization
DTM = maturity-d;
tau = (maturity-d);

if d == maturity    
    % DTM = 0 
    %     fprintf('calc_payoff.m: d==maturity! Yeh!\n');          
    Shat = 1/index_price;
    payoff =  exp(-r*tau) * Shat * max( omega *(index_price - strike), 0);    
else    
    Shat = 1./S;
    payoff =exp(-r*tau) * Shat(:, DTM) .* max(  omega * ( S(:, DTM) - strike), 0);   
end