close all
clear all
clc

%% Part 2
% Battery specifications
% Rhino 6200mAh 4S 50C Lipo Battery Pack w/XT90
% Turnigy Rapid 5500mAh 4S2P 140C Hardcase Lipo Battery Pack w/XT90 Connector

Vmin = 20.8/3.6; %Required speed
Vmax = 15;
V = [5.78:1.31:15]; %Ranges of velocities from minimum required to maximum allowed
voltage = [14.8 14.8];
ah = [6200 5500]*10^(-3);
mass = [705 543]*10^(-3); %kg
totalweight = 6 + mass; %Total weight of UAV
ae = 10; %Aerodynamic efficiency
pe = 0.8; %Propulsion efficiency
g = 9.81;

for j = 1:length(V)
    
    for i = 1:length(ah)
        Espec(i,j) = voltage(i)*ah(i)/mass(i);
        MFbatt(i,j) = mass(i)/totalweight(i);
        E(i,j) = Espec(i,j)*MFbatt(i,j)*pe*ae/(g*V(j))*60^2; %Endurance
        R(i,j) = V(j)*E(i,j); %Range
    end
    
end
figure(1)
plot(V,E(1,:),V,E(2,:))
title('Endurance as a function of cruise velocities');
xlabel('Cruise velocity [m/s]');
ylabel('Endurance [s]');
legend('Rhino 6200mAh', 'Turnigy Rapid 5500mAh')
%% Task 5
%Airfoil NACA 0006, Re = 6*10^6

c_l = [0.22 0.45 0.65 0.83]; %Section lift coefficient
c_d = [0.0055 0.0065 0.0025 0.01]; %Section drag coefficient
ae = [c_l./c_d];

for j = 1:length(ae)
    
    for i = 1:length(ah)
        Emax(i,j) = Espec(i,j)*MFbatt(i,j)*pe*ae(j)/(g*Vmax)*60^2; %Endurance
        Rmax(i,j) = Vmax*Emax(i,j); %Range
    end
    
end

% calculates 
for j = 1:length(ae)
    
    for i = 1:length(ah)
        Emin(i,j) = Espec(i,j)*MFbatt(i,j)*pe*ae(j)/(g*Vmin)*60^2; %Endurance
        Rmin(i,j) = Vmin*Emin(i,j); %Range
    end
    
end

figure(2)
subplot(2,1,1)
plot(ae,Emax(1,:),ae,Emax(2,:))
title('Endurance at Vmax as a function of aerodynamic efficiency');
xlabel('Aerodynamic Efficiency [-]');
ylabel('Endurance [s]');
legend('Rhino 6200mAh', 'Turnigy Rapid 5500mAh')

subplot(2,1,2)
plot(ae,Rmax(1,:),ae,Rmax(2,:))
title('Range at Vmax as a function of aerodynamic efficiency');
xlabel('Aerodynamic Efficiency [-]');
ylabel('Range [m]');

figure(3)
subplot(2,1,1)
plot(ae,Emin(1,:),ae,Emin(2,:))
title('Endurance at Vmin as a function of aerodynamic efficiency');
xlabel('Aerodynamic Efficiency [-]');
ylabel('Endurance [s]');
legend('Rhino 6200mAh', 'Turnigy Rapid 5500mAh')

subplot(2,1,2)
plot(ae,Rmin(1,:),ae,Rmin(2,:))
title('Range at Vmin as a function of aerodynamic efficiency');
xlabel('Aerodynamic Efficiency [-]');
ylabel('Range [m]');

