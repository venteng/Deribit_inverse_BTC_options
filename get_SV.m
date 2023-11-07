% filename: get_SV.m
% written by Huei-Wen Teng: 2021/10/20
function S = get_SV(Z, param, index_price)
% input
% Z = [Z1; Z2]: for vatility, for return  ((3n)x nDay_max)
% param: parameters for the SVCJ model
% index_price: initial stock price
% output
% S

% retrieve n, nDay_max from the data
[nRow, nDay_max] = size(Z);
n = nRow/2;
Z1 = Z(1:n, :);
Z2 = Z((n+1):(2*n), :);

% get the parameters
mu = param(1);
rho = param(2);
alpha = param(3);
beta = param(4);
V0 = param(5);
sig_v = param(6);

% Need to generate all these row-by-row
% 1-2/4
x1 = Z1; % noise for the return yt
x2 = rho* Z1 + sqrt(1-rho^2)*Z2; % noise for the volatiltiy Vt
V = zeros(n, nDay_max);
t = 1;
V(:, 1) = alpha +  beta * V0 + sig_v * sqrt(V0) .* x2(:, t);

for t = 2: nDay_max
    index = find( V(:, (t-1)) < 0 );
    if length(index) > 0
        V(index, (t-1)) = 0;
    end    
    V(:, t) = alpha +  beta * V(:, (t-1)) + sig_v * sqrt(V(:, t-1)) .* x2(:, t);
end

% for t = 1: nDay_max
%     fprintf('t = %g %g \n', t, isreal(V(:, t)));    
%     V(:, t);
% end
Y = zeros(n, nDay_max);
t = 1;
Y(:, t) = mu + sqrt(V0) * x1(:, t);
for t = 2: nDay_max
    Y(:, t) = mu + sqrt(V(:, t-1)) .* x1(:, t);
end
S = index_price * exp(cumsum(Y, 2));