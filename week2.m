clear all
clc

CD0 = 0.035;
CL0 = 0.7;
CL_alpha = 5.7 ;


M = 110 ; % mass [kg]
g = 9.80665 ; % [m/s^2]
W = M*g ; % [N]
S = 2.41 ; % area
K = 0.07 ;

rho = 1.225; % density[kg/m^3]
sigma = 1;

gamma_deg = 0:0.2:5 ; % climb angle
gamma_rad = gamma_deg*pi/180 ;

V = 27:1:35; % velocity [m/s]

%syms T alpha
%f1(T,alpha) = T*cos(alpha)-W*sin(gamma_rad)-0.5*rho*V*V*S*(CD0+K*(CL0 + CL_alpha*alpha)^2);
%f2(T,alpha) = -T*sin(alpha)+W*cos(gamma_rad)-0.5*rho*V*V*S*(CL0+CL_alpha*alpha);

%f1dx(T,alpha) = cos(alpha); %derivative of f wrt x
%f1dy(T,alpha) = -K*S*CL_alpha*rho*(CL0+alpha*CL_alpha)*V^2-T*sin(alpha); %derivative of f wrt y
%f2dx(T,alpha) = -sin(alpha); %derivative of g wrt x
%f2dy(T,alpha) = - T*cos(alpha) - (S*V^2*CL_alpha*rho)/2; %derivative of g wrt y

%x0=200; %initial approximation for x
%y0=0   ; %initial approximation for y
%xprev=0;
%yprev=2;

for j = 1:1:length(gamma_rad)
    for i = 1:1:length(V)

        y0 = 0   ; %initial approximation for y
        x0 = gamma_rad(j)*W + 0.5*rho*sigma*CD0*S*V(i)*V(i) + 2*K*W*W/(rho*sigma*S*V(i)*V(i));  %initial approximation for x

        xprev=0;
        yprev=2;  

        while abs(x0-xprev)>0.03 && abs(y0-yprev)>(1.75*1e-4)

            [F, jacob] = computeF(x0, y0, W, gamma_rad(j), V(i), S,  rho);

            h=-1*inv(jacob)*F;
            xprev=x0;
            yprev=y0;
            x0=x0+h(1);
            y0=y0+h(2);
        end

        t(i,j) = x0;
        a(i,j) = y0;
        p(i,j) = t(i,j)*V(i);
        roc(i,j) = V(i)*sin(gamma_rad(j));

        V_rec(i,j) = V(i);
        gamma_rec(i,j) = gamma_rad(j);


    end
end

% gamma_rec = gamma_rec*57.3;
a = a*57.3;

figure(2);clf
sp(1) = subplot(2,2,1);
surf(gamma_rec, V_rec,t)
xlabel('gamma [rad]')
ylabel('V [m/s]')
zlabel('T [N]')
title('thrust')
grid on ;
hold on

sp(2) = subplot(2,2,2);
surf(gamma_rec, V_rec,a)
grid on
xlabel('gamma [rad]')
ylabel('V [m/s]')
zlabel('alpha [degree]')
title('alpha')

sp(3) = subplot(2,2,3);
surf(gamma_rec, V_rec,p)
grid on
xlabel('gamma [rad]')
ylabel('V [m/s]')
zlabel('P [Nm/s]')
title('power')
zlim([0 15000])

sp(4) = subplot(2,2,4);
surf(gamma_rec, V_rec, roc)
xlabel('gamma [rad]')
ylabel('V [m/s]')
zlabel('RoC')
title('rate of climb')

linkaxes(sp,'x')


function [F, jacob] = computeF(T, alpha, W, gamma_rad, V, S, rho)

f1 = T*cos(alpha)-W*sin(gamma_rad)-0.5*rho*V*V*S*(0.0549 + 0.1326*alpha +0.79363*alpha*alpha +1.6799*alpha*alpha*alpha);
f2 = -T*sin(alpha)+W*cos(gamma_rad)-0.5*rho*V*V*S*(0.6458 + 6.033*alpha -0.54*alpha*alpha -20.966*alpha*alpha*alpha);
f1dx = cos(alpha); %derivative of f1 wrt x
f1dy = - T*sin(alpha) - (S*V^2*rho*((50397*alpha^2)/10000 + (79363*alpha)/50000 + 663/5000))/2; %derivative of f1 wrt y
f2dx = -sin(alpha); %derivative of f2 wrt x
f2dy = (S*rho*((31449*alpha^2)/500 + (27*alpha)/25 - 6033/1000)*V^2)/2 - T*cos(alpha); %derivative of f2 wrt y

% f1 = T*cos(alpha)-W*sin(gamma_rad)-0.5*rho*V*V*S*(0.054916855179320 + 0.132642497237921*alpha +0.935897239571628*alpha*alpha +1.681965269919404*alpha*alpha*alpha);
% f2 = -T*sin(alpha)+W*cos(gamma_rad)-0.5*rho*V*V*S*(0.645801736118242 + 6.034621588677401*alpha -0.556024164979603*alpha*alpha - 20.927471438226657*alpha*alpha*alpha);
% 
% f1dx = cos(alpha); %derivative of f1 wrt x
% f1dy = - T*sin(alpha) - (S*V^2*rho*((2840586811072179*alpha^2)/562949953421312 + (8429812918783719*alpha)/4503599627370496 + 597368701134193/4503599627370496))/2; %derivative of f1 wrt y
% f2dx = -sin(alpha); %derivative of f2 wrt x
% f2dy = (S*rho*((17671678607063301*alpha^2)/281474976710656 + (5008220444422263*alpha)/4503599627370496 - 6794379884522373/1125899906842624)*V^2)/2 - T*cos(alpha); %derivative of f2 wrt y

F = [f1; f2];
jacob = [f1dx f1dy; f2dx f2dy]; %jacobian

end





