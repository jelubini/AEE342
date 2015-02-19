% Joel Lubinitsky
% Fucking around with matlab

clear all
close all
clc

alpha = 0;
n_panels = 8;
m = 2;
p = 4;
tt = 12;

m = m * 0.01;
p = p * 0.1;
tt = tt * 0.01;

zeta = linspace(pi, 2 * pi, n_panels);
xCamber = 0.5 * (1 + cos(zeta));
indexP = max(find(xCamber < p));

xCamber1 = xCamber(1 : indexP);
xCamber2 = xCamber(indexP + 1 : end);
yCamber1 = (m / p ^ 2) * (2 * p * xCamber1 - xCamber1 .^ 2);
yCamber2 = (m / (1 - p) ^ 2) * (1 - (2 * p) + (2 * p * xCamber2) - xCamber2 .^ 2)
yCamber = [yCamber1, yCamber2];

yThickness = (tt / 0.20) * ((0.2969 * sqrt(xCamber)) - (0.1260 .* xCamber) - (0.3516 .* xCamber .^ 2) + (0.2843 .* xCamber .^ 3) - (0.1015 .* xCamber .^ 4));
dydxCamber1 = (2 * m / p) * (1 - xCamber1 / p);
dydxCamber2 = (2 * m / (1 - p) ^ 2) * (p - xCamber2);
dydxCamber = [dydxCamber1, dydxCamber2];
theta = atan(dydxCamber);

xUpper = xCamber - yThickness .* sin(theta);
xLower = xCamber + yThickness .* sin(theta);

yUpper = yCamber + yThickness .* cos(theta);
yLower = yCamber - yThickness .* cos(theta);

X = [xUpper, xLower(end : -1 : 2)];
Y = [yUpper, yLower(end : -1 : 2)];

figure(1)
hold on
axis([-1 2 -1 1])
plot(xCamber, yCamber, 'color', [1 0 0])
plot(X, Y, 'color', [0 0 1])