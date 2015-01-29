%Joel Lubinitsky
%AEE 342 - HW 1
%01/21/15

clear all
close all
clc

%Given
weight = 3200;
areaWing = 180;
density = 0.002377;

%Definitions
forceLift = weight;
velocity = [70 : 1 : 250];

%Calculations
coefficientLift = forceLift ./ (0.5 .* density .* (velocity .^ 2) .* areaWing);
coefficientDrag = 0.0235 + 0.0540 .* (coefficientLift .^ 2);

forceDrag = coefficientDrag .* (0.5 .* density .* (velocity .^ 2) .* areaWing);

ratioLiftDrag = forceLift ./ forceDrag;

ratioLiftDragMax = max(ratioLiftDrag);
velocityMaxLiftDrag = velocity(find(ratioLiftDrag == ratioLiftDragMax))

%Plots
figure(1)
hold on
title('Velocity vs Lift, Drag Coefficient')
xlabel('Velocity [ft/s]')
ylabel('Coefficient')
plot(velocity, coefficientLift, 'color', [1 0 0])
plot(velocity, coefficientDrag, 'color', [0 0 1])
legend('Lift Coefficient', 'Drag Coefficient', 'location', 'best')

figure(2)
hold on
title('Velocity vs Lift/Drag Ratio')
xlabel('Velocity [ft/s]')
ylabel('L/D')
plot(velocity, ratioLiftDrag, 'color', [1 0 0])