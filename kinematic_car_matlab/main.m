clc;clear all;
T = [0:400];
x0 = [0.5;0.5; 0;0];
[t, Y] = ode45(@(t,x) Kinematic_car(t,x), T, x0);

% trajecoty tracking error
figure;
subplot(2,1,1)
plot(t,[Y(:,1)-sin(t/300)])
xlabel('time(s)');ylabel('e_{y1}(m)');
subplot(2,1,2)
plot(t,[Y(:,2)-cos(t/300)])
xlabel('time(s)');ylabel('e_{y2}(m)');

% state 
figure;
subplot(2,2,1)
plot(t,Y(:,1))
xlabel('time (s)');ylabel('y_1 (m)');
subplot(2,2,2)
plot(t,Y(:,2))
xlabel('time (s)');ylabel('y_2 (m)');
subplot(2,2,3)
plot(t,Y(:,3))
xlabel('time (s)');ylabel('\theta (rad)');
subplot(2,2,4)
plot(t,Y(:,4))
xlabel('time (s)');ylabel('\phi (rad)');

% input of system
l=1;
K= [-1 0 0 0 ;0 -1 -2.414 -2.414];
K1 = -1 ;

for i = 1:size(t,1)
theta=Y(i,3);
phi=Y(i,4);
y1=Y(i,1);
y2=Y(i,2);

% desired trajectory 
y1d=sin(t(i)/300);
y2d=cos(t(i)/300);

dy1d=cos(t(i)/300)/300;
ddy1d=-sin(t(i)/300)/90000;
dddy1d=-cos(t(i)/300)/27000000;

dy2d=-sin(t(i)/300)/300;
ddy2d=-cos(t(i)/300)/90000;
dddy2d=sin(t(i)/300)/27000000;


% auxiliary input v1
v1 = dy1d + K1*(y1-y1d);
dv1 = ddy1d + K1*(v1-dy1d);
ddv1 = dddy1d + K1*(dv1-ddy1d);
% dot_y2 = function of (x,v1); Ddot_y2 = (x,v1,dv1)
dy2 = v1*tan(theta);
ddy2 = dv1*tan(theta)+(tan(phi)*v1^2*(tan(theta)^2+1))/(cos(theta)*l);


%u = funcof(x,v1,v2) = D_circ^-1 * (K*e - a_tilde + [dy1d,dddy2d]');
e = [y1-y1d;y2-y2d;dy2-dy2d;ddy2-ddy2d];
av2 = dv1*v1*(tan(phi)*(tan(theta)^2+1))/(l*cos(theta)) + ddv1*tan(theta);% + sin(phi)/(l^2*cos(phi)^2*cos(theta)^5)*(dv1*l*cos(theta)^2*cos(phi)+3*sin(phi)*sin(theta)*v1^2)*v1 ;
a_tilde = [0; av2 ];
d11 = 1/cos(theta);
d12 = 0;
d21 = -(sin(phi)*(dv1*l*cos(theta)^2*cos(phi)+3*sin(phi)*sin(theta)*v1^2))/(l*cos(theta)^2*v1^2);
d22 = cos(theta)^3*l*cos(phi)^2/v1^2;
D_circ_inv = [d11,d12;d21,d22];
u = D_circ_inv*(K*e-a_tilde+[dy1d;dddy2d]);
v(i)=u(1);
dphi(i)=u(2);
end

figure;
subplot(2,1,1)
plot(t,v)
xlabel('time (s)') ;ylabel('v (m/s)')
subplot(2,1,2)
plot(t,dphi)
xlabel('time (s)') ;ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')