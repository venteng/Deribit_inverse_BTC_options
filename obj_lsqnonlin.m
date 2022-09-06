% filename: obj_SVCJ_BRC.m
% updated by Huei-Wen Teng on 2022/1/10
% to calculate the SSE under the BS model using unique option contract
% We use Monte Carlo method to estimate the price
% 1. 
% p is the original prices
% 2. from the unique option contract contract. 
% index_price, cp, d, maturity, strike, ic: ic is from the sorted contract
function f = obj_lsqnonlin(param, model, n, p, U_base, Z_base, index_price_base, omega_base, d_base, maturity_base, strike_base, ic)

% 1. Use iv to calculate ivPrice and see the difference between p and
r = 0; % risk-free rate
y = 0; % dividend rate

m_base = zeros(length(d_base), 1);

for idx = 1:length(d_base)
    m_base(idx,1) = calc_inverse(model, n, U_base, Z_base, omega_base(idx), index_price_base(idx), strike_base(idx), r, y, d_base(idx), maturity_base(idx), param);
end

m = m_base(ic);

f = (p-m);
