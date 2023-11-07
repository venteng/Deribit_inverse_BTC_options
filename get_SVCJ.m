% filename: get_SVCJ.m
% written by Huei-Wen Teng: 2021/10/20
function SVCJ_return  = get_SVCJ(U, Z, param)
% input
% U = [U1;U2]: for jump of size ((2n)x nDay_max)
% Z = [Z1; Z2; Z3]: for jump size in return, for vatility, for return  ((3n)x nDay_max)
% param: parameters for the SVCJ model
% output
% accumulative returns under the SVCJ model

% retrieve n, nDay_max from the data
[nRow, nDay_max] = size(U);
n = nRow/2; % Monte Carlo sample size

U1 = U( 1:n, :);
U2 = U( (n+1) : 2*n, :);
Z1 = Z(1:n, :);
Z2 = Z((n+1):(2*n), :);
Z3 = Z((2*n+1):(3*n), :);

% get the parameters
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


% Need to generate all these row-by-row
% 1/5
J = get_Bernoulli(lambda, U1);
% 2/5
jump2 =  get_exp(1/mu_v, U2); % Jump for the volatiltiy
% 3/5
jump1 = ( mu_y + rho_j * jump2 ) + abs(sig_y) * Z1; % Jump for the return
% 4-5/5
x1 = Z2; % noise for the return yt
x2 = rho* Z2 + sqrt(1-rho^2)*Z3; % noise for the volatiltiy Vt
V = zeros(n, nDay_max);
t = 1;
V(:, 1) = alpha +  beta * V0 + sig_v * sqrt(V0) .* x2(:, t) + jump2(:, t).*J(:, t);

for t = 2: nDay_max
    index = find( V(:, (t-1)) < 0 );
    if length(index) > 0
        V(index, (t-1)) = 0;
    end    
    V(:, t) = alpha +  beta * V(:, (t-1)) + sig_v * sqrt(V(:, t-1)) .* x2(:, t) + jump2(:, t).*J(:, t);
end

% for t = 1: nDay_max
%     fprintf('t = %g %g \n', t, isreal(V(:, t)));    
%     V(:, t);
% end
Y = zeros(n, nDay_max);
t = 1;
Y(:, t) = mu + sqrt(V0) * x1(:, t) + jump1(:,t).*J(:, t);
for t = 2: nDay_max
    Y(:, t) = mu + sqrt(V(:, t-1)) .* x1(:, t) + jump1(:,t).*J(:, t);
end

SVCJ_return = cumsum(Y, 2);