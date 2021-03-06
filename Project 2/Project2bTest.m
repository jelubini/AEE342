% Joel Lubinitsky
% Thickness shit

clear all
close all
clc

n_pan = 10000;
zeta = linspace(pi, 2 * pi, n_pan / 2);
xCamber = 0.5 * (1 + cos(zeta));
tt = 0.15;

yThickness = (tt / 0.20) * ((0.2969 * sqrt(xCamber)) - (0.1260 .* xCamber) - (0.3516 .* xCamber .^ 2) + (0.2843 .* xCamber .^ 3) - (0.1015 .* xCamber .^ 4));

[yThicknessMax, indexTTMax] = max(yThickness)

ttMax = 2 * yThicknessMax

coefficientNew = 0.2969 - 0.1260 - 0.3516 + 0.2843

yThicknessNew = (tt / 0.20) * ((0.2969 * sqrt(xCamber)) - (0.1260 .* xCamber) - (0.3516 .* xCamber .^ 2) + (0.2843 .* xCamber .^ 3) - (coefficientNew .* xCamber .^ 4));

[yThicknessMaxNew, indexTTMaxNew] = max(yThicknessNew)

ttMaxNew = 2 * yThicknessMaxNew


figure(1)
hold on
title('New vs Old Thickness Distribution')
xlabel('X')
ylabel('Y')
axis equal
plot(xCamber, yThickness, '-','color', [0 0 0])
plot(xCamber, -yThickness, '-', 'color', [0 0 0])
plot(xCamber, yThicknessNew, '--', 'color', [0 0 0])
plot(xCamber, -yThicknessNew, '--', 'color', [0 0 0])