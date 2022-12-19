function dx = Kinematic_car(t,x)
y1=x(1);
y2=x(2);
theta=x(3);
phi=x(4);

l=0.1;

% LQR Coefficient
K= [-1 0 0 0 ;0 -1 -2.414 -2.414];
K1 = -1 ;

%trajactory and time derivative of trajectory
y1d=sin(t/300);
y2d=cos(t/300);

dy1d=cos(t/300)/300;
ddy1d=-sin(t/300)/90000;
dddy1d=-cos(t/300)/27000000;

dy2d=-sin(t/300)/300;
ddy2d=-cos(t/300)/90000;
dddy2d=sin(t/300)/27000000;


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
v=u(1);
dphi=u(2);
% v=v1/cos(theta);
% dphi=(1/(v1^2*cos(theta)^2*l)) * v2*l^2*cos(phi)^2*cos(theta)^5 - ddv1*sin(theta)*l^2*cos(phi)^2*cos(theta)^4 - 3*cos(theta)^2*sin(phi)*cos(phi)*l*v1*dv1 + 3*cos(phi)^2*sin(theta)*v1^3 - 3*v1^3*sin(theta);

dx1= v*cos(theta);
dx2= v*sin(theta);
dx3= v*tan(phi)/l;
dx4= dphi;
dx = [dx1;dx2;dx3;dx4];
end

