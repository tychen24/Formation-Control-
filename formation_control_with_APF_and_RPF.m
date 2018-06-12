%% Formation Control + Obstacle Advoidance using APF and RPF 
% Since RPF is nonlinear function, and the graph is unbalanced, the control
% law can't be combined RPF with APF (linear) in an array
global eta d dt iteration 

%plot the directed graph
G = digraph([1 2 2 3 3],[2 1 3 1 4]); 

figure(1)
plot(G)
title('Formation Graph D(G)')

%gain of RPF
eta =10; 

% The distance whcih robots will be in danger of collision (activated zone)
d =6 ; 

%sampling period
dt=0.01;

%numbers of iterations 
iteration =599;

%initialize the array for the state of agents
%z_d = [z1x(0) z1x(k) ...
%       z1y(0) z1y(k) ......
%       .
%       .
%       z8y(0) z8y(k)]
%size(z_d) = [numNodes*2, k]
z_d = zeros(numNodes*2, iteration+1);

%Initialize the array for the distance between each 'adjacent' agents
% dis = [  0   dis12 dis13 dis14;
%        dis21   0   dis23 dis24;
%        dis31 dis32   0   dis34;
%        dis41 dis42 dis43   0  ];
dis = zeros(numNodes,numNodes,length(z_d));

%Initialize the array for the RPF's derivative
% d_v = [  0   dv_12 dv_13 dv_14;
%        dv_21   0   dv_23 dv_24;
%        dv_31 dv_32   0   dv_34;
%        dv_41 dv_42 dv_43   0  ];
d_v =zeros(numNodes,numNodes,length(z_d));

for k =1:iteration+1
    
    z_d(:,1)=z0;
    dis(1,2,k) = norm(z_d(1:2,k) - z_d(3:4,k)); %distance between node 1 and node 2
    dis(1,3,k) = norm(z_d(1:2,k) - z_d(5:6,k)); %distance between node 1 and node 3
    dis(1,4,k) = norm(z_d(1:2,k) - z_d(7:8,k)); %distance between node 1 and node 4
    dis(2,1,k) = dis(1,2,k);                    %distance between node 2 and node 1 (same as dis21)
    dis(2,3,k) = norm(z_d(3:4,k) - z_d(5:6,k)); %distance between node 2 and node 3
    dis(2,4,k) = norm(z_d(3:4,k) - z_d(7:8,k)); %distance between node 2 and node 4
    dis(3,1,k) = dis(1,3,k);                    %distance between node 3 and node 1
    dis(3,2,k) = dis(2,3,k);                    %distance between node 3 and node 2
    dis(3,4,k) = norm(z_d(5:6,k) - z_d(7:8,k)); %distance between node 3 and node 4
    dis(4,1,k) = dis(1,4,k);                    %distance between node 4 and node 1
    dis(4,2,k) = dis(2,4,k);                    %distance between node 4 and node 2
    dis(4,3,k) = dis(3,4,k);                    %distance between node 4 and node 3
    
    d_v(1,2,k) = -4*eta*(d^2 - dis(1,2,k)^2) / (d^2*dis(1,2,k)^5); %potential derivative of node 1 to 2
    d_v(2,1,k) = -d_v(2,1,k);                                      %potential derivative of node 2 to 1
    d_v(1,3,k) = -4*eta*(d^2 - dis(1,3,k)^2) / (d^2*dis(1,3,k)^5); %potential derivative of node 1 to 3
    d_v(3,1,k) = -d_v(1,3,k);
    d_v(1,4,k) = -4*eta*(d^2 - dis(1,4,k)^2) / (d^2*dis(1,4,k)^5);
    d_v(4,1,k) = -d_v(1,4,k);
    d_v(2,3,k) = -4*eta*(d^2 - dis(2,3,k)^2) / (d^2*dis(2,3,k)^5);
    d_v(3,2,k) = -d_v(2,3,k);
    d_v(2,4,k) = -4*eta*(d^2 - dis(2,4,k)^2) / (d^2*dis(2,4,k)^5);
    d_v(4,2,k) = -d_v(2,4,k);
    d_v(3,4,k) = -4*eta*(d^2 - dis(3,4,k)^2) / (d^2*dis(3,4,k)^5);
    d_v(4,3,k) = -d_v(3,4,k);
    
    % if-else condition, make sure that RPF starts to activate when obstacles enter the danger zone
    % if (current distance)^2 - d^2 <=0
    %    de = 1 
    % else
    %    de = 0
    de12 = (dis(1,2,k)^2 - d^2 <= 0); 
    de13 = (dis(1,3,k)^2 - d^2 <= 0);
    de14 = (dis(1,4,k)^2 - d^2 <= 0);
    de23 = (dis(2,3,k)^2 - d^2 <= 0);
    de24 = (dis(2,4,k)^2 - d^2 <= 0);
    de34 = (dis(3,4,k)^2 - d^2 <= 0);
    
    %nonlinear equation of the agent dynamic (z(k+1) = z(k) - APF's*dt - RPF's*dt))
    z_d(1,k+1) = z_d(1,k) - de12*d_v(1,2,k)*dt - de13*d_v(1,3,k)*dt -de14*d_v(1,4,k)*dt...
        - kd*(( (z_d(1,k) - z_d(3,k)) + (z_d(1,k) - z_d(5,k)) ) - c_mat(1,1) )*dt; %z1_x
    
    z_d(2,k+1) = z_d(2,k) - de12*d_v(1,2,k)*dt - de13*d_v(1,3,k)*dt -de14*d_v(1,4,k)*dt...
        - kd*(( (z_d(2,k)-z_d(4,k)) + (z_d(2,k)-z_d(6,k)) ) - c_mat(2,1) )*dt; %z1_y
    
    z_d(3,k+1) = z_d(3,k) - de12*d_v(2,1,k)*dt - de23*d_v(2,3,k)*dt -de24*d_v(2,4,k)*dt...
        - kd*( z_d(3,k) - z_d(1,k) - c_mat(3,1) )*dt; %z2_x
    
    z_d(4,k+1) = z_d(4,k)  - de12*d_v(2,1,k)*dt - de23*d_v(2,3,k)*dt -de24*d_v(2,4,k)*dt...
        - kd*( z_d(4,k) - z_d(2,k) - c_mat(4,1) )*dt; %z2_y
    
    z_d(5,k+1) = z_d(5,k) - de13*d_v(3,1,k)*dt - de23*d_v(3,2,k)*dt -de24*d_v(3,4,k)*dt...
        - kd*( z_d(5,k) - z_d(3,k) - c_mat(5,1) )*dt; %z3_x
    
    z_d(6,k+1) = z_d(6,k) - de13*d_v(3,1,k)*dt - de23*d_v(3,2,k)*dt -de24*d_v(3,4,k)*dt...
        - kd*( z_d(6,k) - z_d(4,k) - c_mat(6,1) )*dt; %z3_y
    
    z_d(7,k+1) = z_d(7,k) - de14*d_v(4,1,k)*dt - de24*d_v(4,2,k)*dt - de34*d_v(4,3,k)*dt ...
        -kd*( z_d(7,k) - z_d(5,k) - c_mat(7,1) )*dt; %z4_x
    
    z_d(8,k+1) = z_d(8,k) - de14*d_v(4,1,k)*dt - de24*d_v(4,2,k)*dt - de34*d_v(4,3,k)*dt ...
        -kd*( z_d(8,k) - z_d(6,k) - c_mat(8,1) )*dt; %z4_y
    
end

figure(2)
plot(z_d(1,1),z_d(2,1),'ro',z_d(3,1),z_d(4,1),'bo',z_d(5,1),z_d(6,1),'go',z_d(7,1),z_d(8,1),'ko')
hold on
for i=2:iteration
    axis([-8 4 -8 4]);
    plot(z_d(1,i),z_d(2,i),'r:.')
    hold on
    plot(z_d(3,i),z_d(4,i),'b:.')
    hold on
    plot(z_d(5,i),z_d(6,i),'g:.')
    hold on
    plot(z_d(7,i),z_d(8,i),'k:.')
    drawnow
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    title('Formation Control and Collision Advoidance using APF and RPF (Digraph) (o=start, x=end)')
    grid on
end
plot(z_d(1,iteration+1), z_d(2,iteration+1),'xr',z_d(3,iteration+1), z_d(4,iteration+1),'xb',...
    z_d(5,iteration+1), z_d(6,iteration+1),'xg',z_d(7,iteration+1), z_d(8,iteration+1),'xk');
legend('Agent1','Agent2','Agent3','Agent4')