% Joel Lubinitsky
% AEE 342 - HW 4.1
% 02/11/15
% Solving for Stagnation Points

clear all
close all
clc

syms x y a

%y == (a ^ 2 - 3 * a * x) / (x ^ 3 + a ^ 2 * x);

x = solve(y == (a ^ 2 - 3 * a * x) / (x ^ 3 + a ^ 2 * x), x)