%Joel Lubinitsky
%AEE 342 - Project 1a: Analysis of Symmetric Airfoils
%RMS Error
%01/30/15

function errorRMS = p1bErrorRMS(nSinks)
nSinks = 99;
ratioPositionChord = [0 0.005 0.0125 0.025 0.050 0.075 0.10 0.20 0.25 0.30 : 0.1 : 0.90 0.95 1.00]';
coefficientPressureExp = [1.000 0.454 0.067 -0.237 -0.450 -0.498 -0.520 -0.510 -0.484 -0.450 -0.369 -0.279 -0.206 -0.132 -0.049 0.055 0.128 1.000]';
t = 0.15;
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

T = 10;
dt = 0.01;
N = (T / dt) + 1;
xyAirfoil = zeros(N, 2);
for i = [-0.001, 0.001]
    xyAirfoil(1, :) = [0, i];
    
    for n = [1 : N - 1]
       xyAirfoil(n + 1, :) = p1bEuler(xyAirfoil(n, :), s, xSink, dt);
    end
end

uAirfoil = 1;
vAirfoil = 0;
for i = [1 : nSinks]
    uAirfoil = uAirfoil + s(i) .* (xyAirfoil(:, 1) - xSink(i)) ./ ((xyAirfoil(:, 1)- xSink(i)) .^ 2 + xyAirfoil(:, 2).^ 2);
    vAirfoil = vAirfoil + s(i) .* xyAirfoil(:, 2)./ ((xyAirfoil(:, 1) - xSink(i)) .^ 2 + xyAirfoil(:, 2).^ 2);
end

velocityFreestream = 0.9996;
[minEndAirfoil, indexEndAirfoil] = min(abs(xyAirfoil(:, 1) - 1));
qAirfoil = sqrt(uAirfoil .^ 2 + vAirfoil .^ 2);
coefficientPressureSim = 1 - (qAirfoil .^ 2) ./ (velocityFreestream .^ 2);

%Root Mean Square Error
sumDiffSq = 0;
for n = [1 : indexEndAirfoil - 1]
    sumDiffSq = sumDiffSq + (coefficientPressureSim(n)  - interp1(ratioPositionChord, coefficientPressureExp, xyAirfoil(n, 1))) .^ 2;
end
errorRMS = sqrt(sumDiffSq ./ (indexEndAirfoil - 1));