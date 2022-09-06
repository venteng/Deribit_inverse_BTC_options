% filename: obj_fminsearch.m
% written by Huei-Wen Teng
% reviewed on 2022/8/31

function SSE = obj_fminsearch(param, model, n, U_base, Z_base, p, index_price_base,  omega_base, d_base, maturity_base, strike_base, ic)

% 1. Use iv to calculate ivPrice and see the difference between p and
r = 0; % risk-free rate
y = 0; % dividend rate

m_base = zeros(length(d_base), 1);

for idx = 1:length(d_base)
    m_base(idx,1) = calc_inverse(model, n, U_base, Z_base, omega_base(idx), index_price_base(idx), strike_base(idx), r, y, d_base(idx), maturity_base(idx), param);
end

m = m_base(ic);
SSE = (p-m)'*(p-m);