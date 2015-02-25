% Joel Lubinitsky
% AEE 342 - HW6: Incompressible Flow over Airfoils (2)
% 02/25/15

clear all
close all
clc

alpha = 4.5;
c     = pi;

m  = 2;
p  = 4;
tt = 12;

alpha = alpha * (pi / 180);
m     = m * 0.01;
p     = p * 0.1;
tt    = tt * 0.01;

theta     = linspace(0, pi, 100);
xCamber  = 0.5 * (1 - cos(theta));
indexP   = find(xCamber < p, 1, 'last');
xCamber1 = xCamber(1 : indexP);
xCamber2 = xCamber(indexP + 1 : end);
zCamber1 = 0.250 .* (0.800 .* xCamber1 - xCamber1 .^ 2);
zCamber2 = 0.111 .* (0.200 + 0.800 .* xCamber2 - xCamber2 .^ 2);
zCamber  = [zCamber1, zCamber2];

% dz1dx = 1 ./ 5 - xCamber1 ./ 2;
% dz2dx = 111 ./ 1250 - (111 .* xCamber2) ./ 500;
% zeta = acos(1 - (2 * p) / c)
% theta = acos(1 - 2 * p)
% A0 = (636491886453715149*pi)/90071992547409920000 + alpha - (797048398351949601*sin(6167402294989009/4503599627370496))/18014398509481984000 + 4915718121213127827985448016935409/405648192073033408478945025720320000
% An = (5734161139222659.*(111.*sin((6167402294989009.*n)./4503599627370496) - 111.*sin(pi.*n)))./(45035996273704960000.*n) - (636491886453715149.*sin((6167402294989009.*n)./4503599627370496 + 6167402294989009./4503599627370496))./(9007199254740992.*(2000.*n + 2000)) - (636491886453715149.*sin((6167402294989009.*n)./4503599627370496 - 6167402294989009./4503599627370496))./(9007199254740992.*(2000.*n - 2000)) - (636491886453715149.*sin(pi.*n))./(9007199254740992.*(2000.*n - 2000)) - (636491886453715149.*sin(pi.*n))./(9007199254740992.*(2000.*n + 2000)) - (5734161139222659.*(sin((6167402294989009.*n)./4503599627370496) - 2.*6.^(1./2).*n.*cos((6167402294989009.*n)./4503599627370496)))./(9007199254740992.*(- 20.*n.^3 + 20.*n))

% New values please fucking work
% A0 = alpha + (35364878569878617467521208754931.*xCamber)./162259276829213363391578010288128 + (636491886453715149.*(5.*xCamber - 2).*(4503599627370496.*pi - 6167402294989009))./202824096036516704239472512860160000 - 35364878569878617467521208754931./405648192073033408478945025720320;
% 
% A = zeros(1, 2);
% for n = [1 2]
% A(n) = (636491886453715149.*(5.*xCamber - 2).*(sin((6167402294989009.*n)./4503599627370496) - sin(pi.*n)))./(22517998136852480000.*n) - (5734161139222659.*sin((6167402294989009.*n)./4503599627370496).*(5.*xCamber - 2))./(90071992547409920.*n);
% end
% int1 = trapz(

% A0 = 

% A0 = alpha + 0.0278 * (acos(0.2) / pi) + 0.0222;

dzdx1 = (1/5) - (1/4) .* (1 - cos(theta(1 : indexP)));
dzdx2 = 0.0888 - 0.111 .* (1 - cos(theta(indexP + 1 : end)));
dzdx = [dzdx1, dzdx2];
dzdxcos1 = dzdx .* cos(theta);
dzdxcos2 = dzdx .* cos(2 * theta);


A0 = alpha - (1/pi) * trapz(theta, dzdx);
A1 = (2/pi) * trapz(theta, dzdxcos1);
A2 = (2/pi) * trapz(theta, dzdxcos2);

gamma0 = 2 .* (A0 .* ((1 + cos(theta)) ./ sin(theta)));
gamma1 = 2 .* (A0 .* ((1 + cos(theta)) ./ sin(theta)) + (A1 .* sin(theta)));
gamma2 = 2 .* (A0 .* ((1 + cos(theta)) ./ sin(theta)) + (A1 .* sin(theta) + A2 .* sin(2 .* theta)));


figure(100)
hold on
% axis([-1 2 -1 1])
plot(theta, zCamber)
plot(theta, dzdx, 'color', 'r')

figure(1)
hold on
axis([0 1 0 2])
plot(xCamber, gamma0, 'color', [1 0 0])
plot(xCamber, gamma1, 'color', [0 1 0])
plot(xCamber, gamma2, 'color', [0 0 1])
legend('\gamma(x) -> A_0', '\gamma(x) -> A_0, A_1', '\gamma(x) -> A_0, A_1, A_2')