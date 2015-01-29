%Joel Lubinitsky
%AEE 342 - Project 1a: Analysis of Symmetric Airfoils
%dy ODE Integrator
%01/23/15

function uv = dxP1a(x, y)

s1 = 0.10;
s2 = -0.07;
s3 = -0.03;

xMin = -4;
xMax = 4;
yMin = -3;
yMax = 3;


uv = zeros(1, 2);
uv(1) = s1 .* (x + 1) ./ ((x + 1) .^ 2 + y .^ 2) + s2 .* x ./ (x .^ 2 + y .^ 2) + s3 .* (x - 1) ./ ((x - 1) .^ 2 + y .^ 2) + 1;
uv(2) = s1 .* y ./ ((x + 1) .^ 2 + y .^ 2) + s2 .* y ./ (x .^ 2 + y .^ 2) + s3 .* y ./ ((x - 1) .^ 2 + y .^ 2);

end