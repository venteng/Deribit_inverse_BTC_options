% filename: get_exp.m
% Written by Huei-Wen Teng, 2021/10/20
% Input
% p: parameter, probability of getting 1
% U: uniform realizations
% Output
% x: ~ Bernoulli(p)
function x = get_exp(lambda, U)

x = -log(U)/lambda;
