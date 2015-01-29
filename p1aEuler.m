%Joel Lubinitsky
%AEE 342 - Project 1a: Analysis of Symmetric Airfoils
%Euler Loop Function
%01/23/15

function xyNext = p1aEuler(xy, dt)
uv = zeros(1, 2);
xyNext = zeros(1, 2);

s1 = 0.10;
s2 = -0.07;
s3 = -0.03;

x = xy(1);
y = xy(2);

uv(1) = s1 .* (x + 1) ./ ((x + 1) .^ 2 + y .^ 2) + s2 .* x ./ (x .^ 2 + y .^ 2) + s3 .* (x - 1) ./ ((x - 1) .^ 2 + y .^ 2) + 1;
uv(2) = s1 .* y ./ ((x + 1) .^ 2 + y .^ 2) + s2 .* y ./ (x .^ 2 + y .^ 2) + s3 .* y ./ ((x - 1) .^ 2 + y .^ 2);

xyNext(1) = uv(1) .* dt + x;
xyNext(2) = uv(2) .* dt + y;
end