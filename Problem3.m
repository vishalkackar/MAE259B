% Problem 3: Simulation of the deformation of elastic beams

%% Clearing

clc; clear all; close all;

%% Globals

global EA EI dL nodes
global dt Pvec v
global mass massMat
global maxIter ep

%% Properties (no gravity or viscous drag forces)

nodes = 50;
edges = nodes - 1;

ti = 0;
tf = 1;
dt = 1e-2;
steps = (tf - ti)/dt + 1;

L = 1;
dL = L/edges;

R = 0.013;
r = 0.011;

E = 70 * 10^9;
I = pi/4 * (R^4 - r^4);

EA = E * pi * (R^2 - r^2);
EI = E * I;

rho = 2700;

% pos y direction is up
P = -2000;
d = 0.75;

maxIter = 150;
charForce = EI/L^2;
ep = charForce * 1e-7;

%% create position vector

q0 = zeros(2*nodes, 1);
for i=1:1:nodes
    q0(2*i-1) = dL * (i-1);
end

%% velocity vector
v = zeros(2*nodes,1);

%% create mass vector and matrix

%equation given in the text
mass = zeros(2*nodes, 1) + pi*(R^2 - r^2)*L*rho/edges;
mass(1) = mass(1)/2;
mass(2) = mass(2)/2;
mass(2*nodes-1) = mass(2*nodes-1)/2;
mass(2*nodes) = mass(2*nodes)/2;
massMat = diag(mass);

%% Create applied force vector

%right now, say that P is at the first point
diff = abs(d - q0(1));
point = 1;

for i=3:2:2*nodes
    if ( abs(d-q0(i)) < diff )
        diff = abs(d - q0(i));
        point = i;
    end
end

Pvec = zeros(2*nodes, 1);
Pvec(point+1) = P;

%% max vertical displacement

maxDisp = zeros(steps,1);

%% iterate

%     q = update4(q0)
%     v = (q-q0)/dt

for i=2:1:steps
    q = update4(q0);
    v = (q-q0)/dt;
    
    q0 = q;
    
    maxDisp(i) = min(q(2:2:length(q)));
    
    figure(1)
    clf();
    x = q0(1:2:length(q0));
    y = q0(2:2:length(q0));
    
    plot(x,y,'ko-')
    xlabel('x [m]')
    ylabel('y [m]')
    title(num2str((i*dt-dt), 'Time = %f'))
    axis equal
    box on
    drawnow
end

%% Plotting

figure(2)
plot(0:dt:tf, maxDisp)
title("Maximum y displacement as a function of time")
xlabel("Time [s]")
ylabel("Displacement [m]")

%% compare against Euler-Bernoulli

c = min(d, L-d);
ym = P * c * (L^2-c^2)^1.5 / (9 * sqrt(3) * EI * L);

fprintf("Predicted max displacement from Euler beam theory: %f m\n", ym)
fprintf("Max displacement from simulation: %f m\n", maxDisp(length(maxDisp)))