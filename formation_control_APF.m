% Formation Control and Obstacle Advoidance strategy 
% wih APF and RPF, Ting-Yang (Gordon) Chen, University of Washington, June 2018

clear all; close all; clc;
global numNodes L L_kro kd z0 iteration dt tspan 

numNodes = 4; %number of agents

%plot the directed graph
G = digraph([1 2 2 3 3],[2 1 3 1 4]); 

figure(1)
plot(G)
title('Formation Graph D(G)')

% Laplacian Matrix of the directed graph
L = [2 -1 -1 0;
    -1 1 0 0;
    0 -1 1 0;
    0 0 -1 1];

% Kronecker product of the laplacian matrix and identity matrix 
L_kro = kron(L,eye(2));

% gain value
kd = 1;

%initial position (z0(1)=x1, z0(2)=y1, z0(3)=x2, z0(4)=y2, ...)
z0 = [ 0 -1, -1 0, -4 -1, 1 -3]'; 

%desired relative position
c21 = [-5 0];
c12 = [5 0];
c31 = [-5 5];
c23 = [0 -5];
c34 = [-5 0];

%desired vectors of relative positions (c_mat = [sum of the desired vectors of positions for x1, y1, x2, y2, x3, y3, x4, y4])
c_mat = [c21(1,1)+c31(1,1), c21(1,2)+c31(1,2), c12(1,1), c12(1,2), c23(1,1), c23(1,2), c34(1,1), c34(1,2)]';

%numbers of iteration (tspan)
iteration = 199;
%% Formation Control using Attractive Potential Function (APF) (ODE45)

interval = 0.01; %time interval
tspan = 0:interval: iteration*interval+2; %total time elapsed

%Control Law: x_dot = -k*(kron(L,eye(2))-c)
%By Kronecker product, we extend the state z to [x1 y1 x2 y2 x3 y3 x4 y4]'
[t,z_ca] = ode45(@(t,z) -kd*((L_kro*z)-c_mat), tspan, z0);

%plot the graph 
figure(2)
plot(z_ca(1,1),z_ca(1,2),'ro',z_ca(1,3),z_ca(1,4),'bo',z_ca(1,5),z_ca(1,6),'go',z_ca(1,7),z_ca(1,8),'ko') %plot the initial node
hold on
for i = 2:iteration
    axis([-8 4 -8 4]);
    plot(z_ca(i,1),z_ca(i,2),'r:.') %plot(x1,y1)
    hold on
    plot(z_ca(i,3),z_ca(i,4),'b:.') %plot(x2,y2)
    hold on
    plot(z_ca(i,5),z_ca(i,6),'g:.') %plot(x3,y3)
    hold on
    plot(z_ca(i,7),z_ca(i,8),'k:.') %plot(x4,y4)
    hold on
    title('Formation Control with APF using ODE45 (o=start, x=end)')
    xlabel('x');
    ylabel('y');
    grid on
    drawnow
end

%plot the node in tfinal in different symbol
plot(z_ca(iteration+1,1), z_ca(iteration+1,2),'xr',z_ca(iteration+1,3), z_ca(iteration+1,4),'xb',...
    z_ca(iteration+1,5), z_ca(iteration+1,6),'xg',z_ca(iteration+1,7), z_ca(iteration+1,8),'xk');
legend('Agent1','Agent2','Agent3','Agent4')

%% Formation Control Using APF (Discrete Differential Equation)

%sampling period
dt=0.01; 
 
%initialize array for the state fo agents
z_da = zeros(numNodes*2, iteration+1);
z_da(:,1)= z0;

%Solve for the first-order differential equation
for k = 1:iteration + 1
    z_da(:,k+1) = z_da(:,k) -(kd*(L_kro*z_da(:,k)-c_mat))*dt;
end

figure(3)
plot(z_da(1,1),z_da(2,1),'ro',z_da(3,1),z_da(4,1),'bo',z_da(5,1),z_da(6,1),'go',z_da(7,1),z_da(8,1),'ko')
hold on
for i=2:iteration
    axis([-8 4 -8 4]);
    plot(z_da(1,i),z_da(2,i),'r:.')
    hold on
    plot(z_da(3,i),z_da(4,i),'b:.')
    hold on
    plot(z_da(5,i),z_da(6,i),'g:.')
    hold on
    plot(z_da(7,i),z_da(8,i),'k:.')
    hold on
    title('Formation Control with APF using discrete differential equation (o=start, x=end)')
    xlabel('x')
    ylabel('y')
    grid on
    drawnow
end
plot(z_da(1,iteration+1), z_da(2,iteration+1),'xr',z_da(3,iteration+1), z_da(4,iteration+1),'xb',...
    z_da(5,iteration+1), z_da(6,iteration+1),'xg',z_da(7,iteration+1), z_da(8,iteration+1),'xk');
legend('Agent1','Agent2','Agent3','Agent4')