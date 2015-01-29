%Joel Lubinitsky
%AEE 342 - Project 1a: Analysis of Symmetric Airfoils
%Euler Loop Function
%01/30/15

function xyNext = p1bEuler(xy, s, xSink, dt)
xyNext = zeros(1, 2);

x = xy(1);
y = xy(2);

u = 1;
v = 0;
for i = [1 : length(s)]
    u = u + s(i) .* (x - xSink(i)) ./ ((x - xSink(i)) .^ 2 + y .^ 2);
    v = v + s(i) .* y ./ ((x - xSink(i)) .^ 2 + y .^ 2);
end


xyNext(1) = u .* dt + x;
xyNext(2) = v .* dt + y;
end