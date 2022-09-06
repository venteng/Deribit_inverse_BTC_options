% filename: calc_delta.m

function delta = calc_delta(model, n, omega, index_price, strike, r, y, d, maturity, param)

switch model
    
    case 'BS'
        
        delta = calc_BS_delta(omega, index_price, strike, r, y, d, maturity, param);
        
    case {'SV', 'SVJ', 'SVCJ'}
        
        % 1. generate common random numbers
        [U, Z] = generate_common_random_numbers(n, (maturity-d));
        % 2. iterate stock prices
        S = iterate_underlying(model, param, index_price, U, Z);
        
        % 3. calculate payoff function
        delta = mean(calc_payoff_delta(S, d, maturity, omega, strike,  r, index_price));
            
    otherwise
        
end
        
        