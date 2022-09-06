% filename: get_Bernoulli.m
% Written by Huei-Wen Teng, 2021/10/20
% Input
% p: parameter, probability of getting 1
% U: uniform realizations
% Output
% x: ~ Bernoulli(p)
function x = get_Bernoulli(p, U)

x = 1.*( U > (1-p));
