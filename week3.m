clc;clear;

M0 = 108; % mass [kg]
rho_fuel = 0.73;    % fuel density [kg/lt]
m_fuel = 22;    % fuel mass [kg]
bingo_fuel_volume = 1;
bingo_fuel = bingo_fuel_volume*rho_fuel;
g = 9.80665 ; % [m/s^2]
W0 = M0*g ; % initial weight [N]
S = 2.41 ; % area

rho0 = 1.225; % density[kg/m^3]
T0 = 288.16; % temperature[K]
B = 0.00650; % lapse rate[K/m]
R = 287; % gas constant[J/kgK]
sigma = 1;

% gamma_deg = 5 ; % climb angle
% gamma_rad = gamma_deg*pi/180 ;

D = 0.53; % diameter of the propeller [m]
V_IAS = 27; % velocity (indicated air speed) [m/s]
h0 = 0; % initial altitude [m]
hmax = 18000*0.3048; %[m]
dt = 60; % [s]

W = W0;
consumed_fuel = 0;
V_TAS = V_IAS;
h = h0;
rho=rho0;
t_climb = 0;
t_cruise = 0;
j=1;

load_data = importdata("load_data_new.mat");
rpm_data = importdata("rpm_data_new.mat");
power_data = importdata("power_data_new.mat");
torque_data = importdata("torque_data_new.mat");
fuelconsumption_data = importdata("fuelconsumption_data_new.mat");
load_rpm_vs_torque = importdata('load_rpm_vs_torque.mat');

[load_mesh,rpm_mesh] = meshgrid(load_data,rpm_data);
F1 = scatteredInterpolant(rpm_mesh(:), power_data(:), load_mesh(:),'linear','none');

power_list = 0:0.2:11;
d_power = power_list(2)-power_list(1);

rpm_list = min(rpm_data):50:max(rpm_data);
d_rpm = rpm_list(2)-rpm_list(1);


% [loadmesh, rpmmesh] = meshgrid(load_data,rpm_data);
% F2 = scatteredInterpolant(loadmesh(:), rpmmesh(:), torque_data(:),'linear','none');


[power_mesh, rpm_mesh2] = meshgrid(power_list, rpm_list);

% [loadmesh, rpmmesh] = meshgrid(load_data,rpm_data);
% F3 = scatteredInterpolant(loadmesh(:), rpmmesh(:), fuelconsumption_data(:),'linear','none');

load_list = 0:0.1:100;
d_load = load_list(2)-load_list(1);

[loadmesh2, rpmmesh2] = meshgrid(load_list, rpm_list);


gamma_rad = 0;
rho = rho0*(1-B*hmax/T0)^(g/(R*B))*(T0/(T0-B*hmax));
V_TAS = V_IAS*sqrt(rho0/rho);

N = 1000;
Range_final = nan(N,1);

tic
while m_fuel > bingo_fuel

    W1 = W;
    W = W - consumed_fuel*g;
    W2 = W;

    y0 = 0   ; %initial approximation for alpha
    x0 = gamma_rad*W + 0.5*rho*sigma*S*V_TAS*V_TAS*0.0549;  % initial approximation for thrust

    xprev=0;
    yprev=2;

    while abs(x0-xprev)>0.03 && abs(y0-yprev)>(1.75*1e-4)

        [F, jacob] = computeF(x0, y0, W, gamma_rad, V_TAS, S,  rho);

        k=-1*(jacob\F);
        xprev=x0;
        yprev=y0;
        x0=x0+k(1);
        y0=y0+k(2);
    end

    t = x0; % thrust
    a = y0*180/pi; % alpha
    p = t*V_TAS; % required power
    pr_SL = p/(1.132*rho/rho0-0.132); % required power aat sea level



    err = 1;
    RPM0 = (V_TAS*60/(D*0.7));
    RPM = 1000;
    while err > 0.1
%         J = V_TAS/((RPM0/60)*D);
%         ct = -0.1543*J*J + 0.0162*J + 0.2125;
%         T = ct*rho*(RPM0/60)^2*D^4;
% 
%         RPM = sqrt(T/(ct*rho*D^4)*3600);
%         err = abs(RPM-RPM0)
%         RPM0 = RPM;

        % RPM = RPM0 - f/df
        RPM = RPM0 - (rho*(RPM0/60)^2*D^4*(-0.1543*(V_TAS/((RPM0/60)*D))*(V_TAS/((RPM0/60)*D)) + 0.0162*(V_TAS/((RPM0/60)*D)) + 0.2125) - t)/((D^4*RPM0*rho*((243*V_TAS)/(250*D*RPM0) - (13887*V_TAS^2)/(25*D^2*RPM0^2) + 17/80))/1800 - (D^4*RPM0^2*rho*((243*V_TAS)/(250*D*RPM0^2) - (27774*V_TAS^2)/(25*D^2*RPM0^3)))/3600);  
        %RPM0 - (rho*(RPM/60)^2*D^4*(-0.1543*(V_TAS/((RPM/60)*D))*(V_TAS/((RPM/60)*D)) + 0.0162*(V_TAS/((RPM/60)*D)) + 0.2125) - t)/((D^4*RPM*rho0*((243*V_TAS)/(250*D*RPM) - (13887*V_TAS^2)/(25*D^2*RPM^2) + 17/80))/1800 - (D^4*RPM^2*rho*((243*V_TAS)/(250*D*RPM^2) - (27774*V_TAS^2)/(25*D^2*RPM^3)))/3600);
        %((rho*(RPM/60)^2*D^4*(-0.1543*(V_TAS/((RPM/60)*D))*(V_TAS/((RPM/60)*D)) + 0.0162*(V_TAS/((RPM/60)*D)) + 0.2125) - t)/((D^3*rho*(972*V_TAS + 425*D*RPM))/3600000)); % RPM0 - ((rho*(RPM/60)^2*D^4*(-0.1543*(V_TAS/((RPM/60)*D))*(V_TAS/((RPM/60)*D)) + 0.0162*(V_TAS/((RPM/60)*D)) + 0.2125) - t)/((D^3*rho0*(972*V_TAS + 425*D*RPM))/3600000));
        err = abs(RPM-RPM0);
        RPM0 = RPM;
     end
    rpm_final = RPM;
    J = V_TAS/((rpm_final/60)*D);
    nu = -8.7193*J^6 + 26.276*J^5 - 30.199*J^4 + 16.955*J^3 - 5.9498*J^2 + 2.5107*J - 0.0133;
    ct = -0.1543*J^2 + 0.0162*J + 0.2125;
    Range  = V_TAS/ct * + 0.6458/0.0549 * log(W1/W2); %V_TAS/(ct*sqrt(K*CD0))*(atan(2/(rho*V_TAS*V_TAS*S)*sqrt(K/CD0)*W1)-atan(2/(rho*V_TAS*V_TAS*S)*sqrt(K/CD0)*W2));
    Range_final(j,1) = Range;
    ps = pr_SL/nu; % shaft power at sea level

    load_mesh2 = F1(rpm_mesh2(:), power_mesh(:));
    load_mesh2 = reshape(load_mesh2, size(power_mesh));

    power = ps;
    x = power/1000;

    x1 = floor(x/d_power)+1;
    x2 = ceil(x/d_power)+1;

    rpm =rpm_final;
    y = rpm;

    y1 = floor((y-min(rpm_list))/d_rpm)+1;
    y2 = ceil((y-min(rpm_list))/d_rpm)+1;
    
    %     bilinear interpolation
    load_final = 1/((power_list(x2) - power_list(x1))*(rpm_list(y2) - rpm_list(y1)))*[(power_list(x2)-x) (x-power_list(x1))] *[load_mesh2(y1 ,x1) load_mesh2(y2,x1); load_mesh2(y1,x2) load_mesh2(y2,x2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];

    %     torquemesh = F2(loadmesh2(:), rpmmesh2(:));
    %     torquemesh = reshape(torquemesh, size(rpmmesh2));
    E = load_rpm_vs_torque.coeff_d;
    C = load_rpm_vs_torque.coeff;
    load_max = load_rpm_vs_torque.load_max;
    rpm_max = load_rpm_vs_torque.rpm_max;

    l = load_list/load_max;
    r = rpm_list/rpm_max;
    [l_mesh, r_mesh] = meshgrid(l, r);

    %     surface fit
    torquemesh = C(1) + C(2)*l_mesh + C(3)*r_mesh + C(4)*l_mesh.^2 + C(5)*r_mesh.^2 + C(6)*l_mesh.*r_mesh;
    fuelconsmesh = E(1) + E(2)*l_mesh + E(3)*r_mesh + E(4)*l_mesh.^2 + E(5)*r_mesh.^2 + E(6)*l_mesh.*r_mesh;
    

    load = load_final;
    z = load;

    z1 = floor(z/d_load)+1;
    z2 = ceil(z/d_load)+1;

    %  bilinear interpolation
    torque_final = 1/((load_list(z2) - load_list(z1))*(rpm_list(y2) - rpm_list(y1)))*[(load_list(z2)-z) (z-load_list(z1))] *[torquemesh(y1 ,z1) torquemesh(y2,z1); torquemesh(y1,z2) torquemesh(y2,z2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];
    BSFC_final = 1/((load_list(z2) - load_list(z1))*(rpm_list(y2) - rpm_list(y1)))*[(load_list(z2)-z) (z-load_list(z1))] *[fuelconsmesh(y1 ,z1) fuelconsmesh(y2,z1); fuelconsmesh(y1,z2) fuelconsmesh(y2,z2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];

    consumed_fuel = BSFC_final/(3600*10^6)*ps*(1.132*rho/rho0-0.132)*dt; % [Newton]

    m_fuel = m_fuel - consumed_fuel;
    t_cruise = t_cruise + dt;
    j=j+1;

end

m_fuel
t_cruise/3600
toc


% surf(loadmesh, rpmmesh, torque_data)

% Newton Raphson
function [F, jacob] = computeF(T, alpha, W, gamma_rad, V, S, rho)

f1 = T*cos(alpha)-W*sin(gamma_rad)-0.5*rho*V*V*S*(0.0549 + 0.1326*alpha +0.79363*alpha*alpha +1.6799*alpha*alpha*alpha);
f2 = -T*sin(alpha)+W*cos(gamma_rad)-0.5*rho*V*V*S*(0.6458 + 6.033*alpha -0.54*alpha*alpha -20.966*alpha*alpha*alpha);
f1dx = cos(alpha); %derivative of f1 wrt x
f1dy = - T*sin(alpha) - (S*V^2*rho*((50397*alpha^2)/10000 + (79363*alpha)/50000 + 663/5000))/2; %derivative of f1 wrt y
f2dx = -sin(alpha); %derivative of f2 wrt x
f2dy = (S*rho*((31449*alpha^2)/500 + (27*alpha)/25 - 6033/1000)*V^2)/2 - T*cos(alpha); %derivative of f2 wrt y


F = [f1; f2];
jacob = [f1dx f1dy; f2dx f2dy]; %jacobian

end