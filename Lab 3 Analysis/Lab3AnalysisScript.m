clc; clear; close all;
%% Top Level Variables
controller_definition_file_name = "Good_Controller.mat";
controller_data_file_name = "jackson_controller.csv";

% Plant Definition
R = .00635;
k=3.7;
Mc=.915;
Mr=0.210;
Imot=3.87e-7;
L=0.6;
g=9.8;

a1= 1.1722; 
a2= -7.9096; 

b11=Mc + Mr + Imot*k^2/R^2;
b12=-.5*Mr*L;
b22=2*L/3;
nump=-a1*[(2*L/3-b22) 0 g];
denp=[b11*b22+b12, -a2*b22, -g*b11, g*a2, 0];

G=tf(nump,denp);

%% Load in Data

% Load in ZPK controller definition
load(controller_definition_file_name);
C = zpk(z, p, k)

% Load in data
controller_data = readtable(controller_data_file_name);

%% Plot and save figures

% Voltage
voltage = figure();
hold on;
plot(controller_data.SimulationTime - controller_data.SimulationTime(1), controller_data.OutputVoltage);
xlabel("Elapsed Time (sec)");
ylabel("Output Voltage (V)");
title("Output Voltage vs Time");
xlim([0 10]);
saveas(voltage, "Voltage_vs_Time_" + strrep(controller_data_file_name, ".csv", "") + ".png");

% Position
position = figure();
hold on;
plot(controller_data.SimulationTime - controller_data.SimulationTime(1), controller_data.CartPosition);
xlabel("Elapsed Time (sec)");
ylabel("Cart Position (cm)");
title("Cart Position vs Time");
xlim([0 10]);
saveas(position, "Position_vs_Time_" + strrep(controller_data_file_name, ".csv", "") + ".png");

% Angle
angle = figure();
hold on;
plot(controller_data.SimulationTime - controller_data.SimulationTime(1), controller_data.PendulumAngle);
xlabel("Elapsed Time (sec)");
ylabel("Pendulum Angle (rad)");
title("Pendulum Angle vs Time");
xlim([0 10]);
saveas(angle, "PendulumAngle_vs_Time_" + strrep(controller_data_file_name, ".csv", "") + ".png");

%% Plot Open Loop and Compensator Figures
L = C*G;
loop_bode = figure();
bode(L);
saveas(loop_bode, "Loop_bode" + strrep(controller_data_file_name, ".csv", "") + ".png");

loop_step = figure();
step(feedback(L,1));
saveas(loop_step, "Loop_step" + strrep(controller_data_file_name, ".csv", "") + ".png");

loop_root_locus = figure();
rlocus(L);
saveas(loop_root_locus, "Loop_root_locus" + strrep(controller_data_file_name, ".csv", "") + ".png");

compensator_bode = figure();
bode(C);
saveas(compensator_bode, "Compensator_bode" + strrep(controller_data_file_name, ".csv", "") + ".png");

%% Plot Plant Root Locus
plant_rl = figure();
rlocus(G);
saveas(plant_rl, "plant_root_locus.png");






