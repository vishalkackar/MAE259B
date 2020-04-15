% Problem 2: N-connected spheres falling inside viscous fluid

%% Clear variables and figures

clc;
clear all;
close all;

%% Globals

global ep maxIter
global EA EI L dL
global nodes middle Rvec
global q0 dt v
global mu Weight mass massMat
global R Rm
%% Define and create necessary properties

%nodes = input("Number of nodes: ");
nodes = 21;
links = nodes - 1;
middle = ceil(nodes/2);

dt = .01;
t = 0;
tf = 50;
steps = (tf - t)/dt;

L = .1;
dL = L/links;

R = dL/10;
Rm = .025;

Rvec = zeros(2*nodes,1);
for i=1:1:nodes
    if (i == middle)
        Rvec(2*i-1) = Rm;
        Rvec(2*i) = Rm;
    else
        Rvec(2*i-1) = R;
        Rvec(2*i) = R;
    end
end

rho_m = 7000;
rho_f = 1000;
rho_diff = (rho_m - rho_f);

g = -9.81;

mu = 1000;

Y = 1 *10^9;
r0 = .001;
EA = Y*pi*r0^2;
EI = Y/4*pi*r0^4;

%% create mass vector and matrix

mass = zeros(2*nodes,1);
for i=1:1:nodes
    mass(2*i-1) = 4/3 * pi * rho_m * Rvec(2*i-1)^3;
    mass(2*i) = 4/3 * pi * rho_m * Rvec(2*i)^3;
end

massMat = diag(mass);

%% create weight vector, no weight in x direction
Weight = zeros(2*nodes, 1);
for i=1:1:nodes
    Weight(2*i) = 4/3 * pi * Rvec(2*i)^3 * rho_diff * g;
end

q0 = zeros(2*nodes,1);

for i=1:1:nodes
    q0(2*i-1) = L*(i-1)/links;
end

v = zeros(2*nodes,1);

midpos = zeros(steps,1);
midvel = zeros(steps,1);

%% Error checking info

Fchar = EI/dL^2;
ep = Fchar * 1e-6;
maxIter = 150;

%% video set up
% videoFile = VideoWriter('21-Node-Simulation.mp4','MPEG-4');
% videoFile.FrameRate = 30;
% open(videoFile);
% k = 0;

%% run the iterations

for i=2:1:steps
    q = update2(q0);
    v = (q-q0)/dt;
    
    q0 = q;
    
    midpos(i) = q(middle*2);
    midvel(i) = v(middle*2);
    
%     if(mod(i,5) == 0)
%         figure(1)
%         clf();
%         x = q(1:2:length(q));
%         y = q(2:2:length(q));
% 
%         plot(x,y,'ro-')
%         xlabel('x [m]')
%         ylabel('y [m]')
%         title(num2str((i*dt), 'Time = %f'))
%         axis equal
%         box on
%         drawnow
        
%         FR(k+1) = getframe(gcf);
%         writeVideo(videoFile, FR(k+1));
%         k = k+1;
%     end
    
end

% close(videoFile);

%% plot middle sphere

figure(2)
plot(linspace(0,tf,steps), midpos)
title("Position of the middle sphere vs time")
xlabel("Time [s]")
ylabel("Position [m]")

figure(3)
plot(linspace(0,tf,steps), midvel)
title("Velocity of the middle sphere vs time")
xlabel("Time [s]")
ylabel("Velocity [m/s]")


x = q(1:2:length(q));
y = q(2:2:length(q));

figure(4)
plot(x,y,'ko-')
axis equal
box on
title("Final shape of the beam")
xlabel("x [m]")
ylabel("y [m]")


%% misc info

fprintf("Terminal velocity of the middle node: %f m/s\n",midvel(length(midvel)));