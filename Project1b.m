%Joel Lubinitsky
%AEE 342 - Project 1b: Analysis of Symmetric Airfoils
%01/30/15

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

%Pressure Coefficients (NACA)


%% Calculations
nSinks = 99;

xSink = zeros(1, nSinks);
for i = [1 : length(xSink)]
    xSink(i) = i / (nSinks + 1);
end

yAirfoil = (t ./ 0.20) .* (0.2969 .* sqrt(xSink) - 0.1260 .* xSink - 0.3516 .* xSink .^ 2 + 0.2843 .* xSink .^ 3 - 0.1015 .* xSink .^ 4);

M = zeros(nSinks, nSinks);
for j = [1 : nSinks]
    for i = [1 : nSinks]
        M(j, i) = atan2(yAirfoil(j), (xSink(j) - xSink(i)));
    end
end

R = -yAirfoil';
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
uv = zeros(N, 2);

%Run Loop
figure(1)
hold on
axis([xMin xMax yMin yMax])

for i = [1 : 20]
    xy(1, :) = [x(1), y(i)];
    uv(1, :) = [u(1), v(i)];
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
for i = [-0.001, 0.001]
    xyAirfoil(1, :) = [0, i];
    
    for n = [1 : N - 1]
       xyAirfoil(n + 1, :) = p1bEuler(xyAirfoil(n, :), s, xSink, dt);
    end
    
    plot(xyAirfoil(:, 1), xyAirfoil(:, 2), 'color', [1 0 0])
end
%% Plots
figure(2)
hold on
axis([xMin xMax yMin yMax])
contourf(x, y, v)