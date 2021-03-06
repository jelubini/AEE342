%Joel Lubinitsky
%AEE 342 - Project 1c: Analysis of Symmetric Airfoils
%02/06/15

clear all
close all
clc

%% Known
%Domain
xMin = -1;
xMax = 2;

yMin = -1;
yMax = 1;

%NACA0015
t = 0.15;

%Pressure Coefficients (NACA0015)
ratioPositionChord = [0 0.005 0.0125 0.025 0.050 0.075 0.10 0.20 0.25 0.30 : 0.1 : 0.90 0.95 1.00]';
coefficientPressureExp = [1.000 0.454 0.067 -0.237 -0.450 -0.498 -0.520 -0.510 -0.484 -0.450 -0.369 -0.279 -0.206 -0.132 -0.049 0.055 0.128 1.000]';

%% Calculations
nSinks = 99;

xSink = zeros(1, nSinks);
for i = [1 : length(xSink)]
    xSink(i) = i / (nSinks + 1);
end

yAirfoil = (t ./ 0.20) .* (0.2969 .* sqrt(xSink) - 0.1260 .* xSink - 0.3516 .* xSink .^ 2 + 0.2843 .* xSink .^ 3 - 0.1015 .* xSink .^ 4);
dydxAirfoil = (t ./ 0.20) .* ((0.14845 .* xSink .^ -0.5) - (0.1260) - (0.7032 .* xSink) + (0.8529 .* xSink .^ 2) - (0.4060 .* xSink .^ 3));

M = zeros(nSinks, nSinks);
for j = [1 : nSinks]
    for i = [1 : nSinks]
        M(j, i) = (yAirfoil(j) / ((xSink(j) - xSink(i)) ^ 2 + yAirfoil(j) ^ 2)) - (((xSink(j) - xSink(i)) / ((xSink(j) - xSink(i)) ^ 2 + yAirfoil(j) ^ 2)) * ((t / 0.20) * ((0.14845 * xSink(j) ^ -0.5) - (0.1260) - (0.7032 * xSink(j)) + (0.8529 * xSink(j) ^ 2) - (0.4060 * xSink(j) ^ 3))));
    end
end

R = dydxAirfoil';
s = M\R;

[x, y] = meshgrid(linspace(xMin, xMax, 30), linspace(yMin, yMax, 20));

u = 1;
v = 0;
for i = [1 : nSinks]
    u = u + s(i) .* (x - xSink(i)) ./ ((x - xSink(i)) .^ 2 + y .^ 2);
    v = v + s(i) .* y ./ ((x - xSink(i)) .^ 2 + y .^ 2);
end

%Initialize Loop
T = 10;
dt = 0.01;
N = (T / dt) + 1;
xy = zeros(N, 2);

%Run Loop
figure(1)
hold on
title('Field Plot displaying Airfoil, Sources, Streamlines [nSources = 99]')
xlabel('X')
ylabel('Y')
axis([xMin xMax yMin yMax])

for i = [1 : 20]
    xy(1, :) = [x(1), y(i)];
    for n = [1 : N - 1]
        xy(n + 1, :) = p1bEuler(xy(n, :), s, xSink, dt);
    end
    plot(xy(:, 1), xy(:, 2))
end

plot(xSink, yAirfoil, '*', 'color', [1 0 0])
plot(xSink, -yAirfoil, '*', 'color', [1 0 0])
plot(xSink, 0, 'o', 'color', [0 1 0])

%Airfoil Streamlines
xyAirfoil = zeros(N, 2);
for i = [-yAirfoil(1), yAirfoil(1)]
    xyAirfoil(1, :) = [xSink(1), i];
    
    for n = [1 : N - 1]
       xyAirfoil(n + 1, :) = p1bEuler(xyAirfoil(n, :), s, xSink, dt);
    end
    
    plot(xyAirfoil(:, 1), xyAirfoil(:, 2), 'color', [1 0 1])
end

%Pressure Coefficients
velocityFreestream = sqrt(mean(u(:, 1)) .^ 2 + mean(v(:, 1)) .^ 2);
[minEndAirfoil, indexEndAirfoil] = min(abs(xyAirfoil(:, 1) - 1));

uAirfoil = 1;
vAirfoil = 0;
for i = [1 : nSinks]
    uAirfoil = uAirfoil + s(i) .* (xyAirfoil(:, 1) - xSink(i)) ./ ((xyAirfoil(:, 1)- xSink(i)) .^ 2 + xyAirfoil(:, 2).^ 2);
    vAirfoil = vAirfoil + s(i) .* xyAirfoil(:, 2)./ ((xyAirfoil(:, 1) - xSink(i)) .^ 2 + xyAirfoil(:, 2).^ 2);
end

qAirfoil = sqrt(uAirfoil .^ 2 + vAirfoil .^ 2);
coefficientPressureSim = 1 - (qAirfoil .^ 2) ./ (velocityFreestream .^ 2);

%Root Mean Square Error
nSinkValues = [3 5 9 19 29 49 99];
errorRMS = zeros(length(nSinkValues), 1);
for n = [1 : length(nSinkValues)]
    errorRMS(n) = p1bErrorRMS(nSinkValues(n));
end

%% Plots
%Source-Sink Distribution
figure(2)
hold on
title('Source Position vs Strength [nSources = 99]')
xlabel('Position (x/c)')
ylabel('Source Strength')
plot(xSink, s, '-*')

%Pressure Coefficients
figure(3)
hold on
axis([0 1 -1 1])
title('Pressure Coefficient Distribution [nSources = 99]')
xlabel('Position (x/c)')
ylabel('Cp')
plot(ratioPositionChord, coefficientPressureExp, 'o', 'color', [1 0 0])
plot(xyAirfoil(:, 1), coefficientPressureSim)

%RMS Error
figure(4)
hold on
title('Root Mean Square Error Convergence')
xlabel('Number of Sources')
ylabel('Error Cp')
loglog(nSinkValues, errorRMS, '-o')