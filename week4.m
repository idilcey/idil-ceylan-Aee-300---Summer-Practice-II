clc;clear;

M0 = 108 ; % mass [kg]
m_fuel = 22;    % fuel mass [kg]
bingo_fuel = 2;
g = 9.80665 ; % [m/s^2]
W0 = M0*g ; % [N]
S = 2.41 ; % area

rho0 = 1.225; % density[kg/m^3]
T0 = 288.16; % temperature[K]
B = 0.00650; % lapse rate[K/m]
R = 287; % gas constant[J/kgK]

sigma = 1;

gamma_deg = 5 ; % climb angle
gamma_rad = gamma_deg*pi/180 ;

D = 0.53; % diameter of the propeller [m]
V_IAS = 27; % velocity [m/s]
h0 = 0; % initial altitude [m]
hmax = 18000*0.3048; %[m]
dt = 60; % [s]


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




[power_mesh, rpm_mesh2] = meshgrid(power_list, rpm_list);



load_list = 0:0.1:100;
d_load = load_list(2)-load_list(1);

[loadmesh2, rpmmesh2] = meshgrid(load_list, rpm_list);





W = W0;
consumed_fuel = 0;
V_TAS = V_IAS;
h = h0;
rho=rho0;
t_climb = 0;
t_cruise = 0;
j=1;


while h < hmax

    [W, h, rho, V_TAS] = myfun(consumed_fuel, g, V_TAS, gamma_rad, rho0, B,  T0, R, V_IAS,W, h, dt);
    
    CD = 0.0549 + 0.1326*gamma_rad + 0.9359*gamma_rad^2 + 1.6820*gamma_rad^3;

    y0 = 0   ; %initial approximation for y
    x0 = gamma_rad*W + 0.5*rho*sigma*S*V_TAS*V_TAS*CD; %CD0 + 2*K*W*W/(rho*sigma*S*V_TAS*V_TAS);  %initial approximation for x

    xprev=0;
    yprev=2;

    while abs(x0-xprev)>0.03 && abs(y0-yprev)>(1.75*1e-4)

        [F, jacob] = computeF(x0, y0, W, gamma_rad, V_TAS, S,  rho);

        k=-1*inv(jacob)*F;
        xprev=x0;
        yprev=y0;
        x0=x0+k(1);
        y0=y0+k(2);
    end

    t = x0;
    a = y0*180/pi;
    p = t*V_TAS;
    pr_SL = p/(1.132*rho/rho0-0.132); %p*sqrt(rho0/rho); required power at sea level
    
    syms RPM fun
    fun = (17*rho*D^4*RPM^2)/288000 + (27*V_TAS*rho*D^3*RPM)/100000 - (1543*V_TAS^2*rho*D^2)/10000-t;
    sol = solve(fun, RPM);
    sol = double(sol);
    for j=1:2
    if (sol(j)>0)
        RPM=sol(j);
    end
    end

    J = V_TAS/((RPM/60)*D);
    ct = -0.1543*J*J + 0.0162*J + 0.2125;
    t_new = (ct*rho*(RPM/60)^2*D^4); 
    cp = -0.149810925027*J^3 + 0.074899811242*J^2 + 0.055465658837*J + 0.107808966045;
    ps_new = cp*rho*(RPM/60)^3*D^5;

    rpm_final = RPM;
    J = V_TAS/((rpm_final/60)*D);
    nu = -8.7193*J^6 + 26.276*J^5 - 30.199*J^4 + 16.955*J^3 - 5.9498*J^2 + 2.5107*J - 0.0133;
    ps = pr_SL/nu;
    ps_alt = ps*(1.132*rho/rho0-0.132);
    %     ps(i) = ps(i)/1000; % [kW]

    load_mesh2 = F1(rpm_mesh2(:), power_mesh(:));
    load_mesh2 = reshape(load_mesh2, size(power_mesh));
% 
%     power_mesh2 = F2(rpm_mesh2(:), loadmesh2(:));
%     power_mesh2 = reshape(power_mesh2, size(loadmesh2));
    


    power =ps;
    x = power/1000;
    for i = 1:length(power_list)-1
        if x > power_list(i)  &&  x< power_list(i+1)
            x1 = i;
            x2 = i+1;
        end
    end



    rpm =rpm_final;
    y = rpm;

    for i = 1:length(rpm_list)-1
        if y > rpm_list(i) && y < rpm_list(i+1)
            y1 = i;
            y2 = i+1;
        end
    end
    
      load_final = 1/((power_list(x2) - power_list(x1))*(rpm_list(y2) - rpm_list(y1)))*[(power_list(x2)-x) (x-power_list(x1))] *[load_mesh2(y1 ,x1) load_mesh2(y2,x1); load_mesh2(y1,x2) load_mesh2(y2,x2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];

    
    E = load_rpm_vs_torque.coeff_d;
    C = load_rpm_vs_torque.coeff;
    load_max = load_rpm_vs_torque.load_max;
    rpm_max = load_rpm_vs_torque.rpm_max;

    l = load_list/load_max;
    r = rpm_list/rpm_max;
    [l_mesh, r_mesh] = meshgrid(l, r);


    torquemesh = C(1) + C(2)*l_mesh + C(3)*r_mesh + C(4)*l_mesh.^2 + C(5)*r_mesh.^2 + C(6)*l_mesh.*r_mesh;
    fuelconsmesh = E(1) + E(2)*l_mesh + E(3)*r_mesh + E(4)*l_mesh.^2 + E(5)*r_mesh.^2 + E(6)*l_mesh.*r_mesh;


    load =load_final;
    z = load;
    for i = 1:length(load_list)-1
        if z > load_list(i) && z < load_list(i+1)
            z1 = i;
            z2 = i+1;
        end
    end

    torque_final = 1/((load_list(z2) - load_list(z1))*(rpm_list(y2) - rpm_list(y1)))*[(load_list(z2)-z) (z-load_list(z1))] *[torquemesh(y1 ,z1) torquemesh(y2,z1); torquemesh(y1,z2) torquemesh(y2,z2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];
    BSFC_final = 1/((load_list(z2) - load_list(z1))*(rpm_list(y2) - rpm_list(y1)))*[(load_list(z2)-z) (z-load_list(z1))] *[fuelconsmesh(y1 ,z1) fuelconsmesh(y2,z1); fuelconsmesh(y1,z2) fuelconsmesh(y2,z2)]*[(rpm_list(y2)-y); (y-rpm_list(y1))];

    consumed_fuel = BSFC_final/(3600*10^6)*ps*(1.132*rho/rho0-0.132)*dt; % [Newton]
    
    if consumed_fuel < 0 || load_final > 100
        break
    end

    m_fuel = m_fuel - consumed_fuel;
    t_climb = t_climb + dt;
    tclimb(j) = t_climb/3600;
    loadfinal(j) = load_final;

    j=j+1;
    
end






function [W, h, rho, V_TAS] = myfun(consumed_fuel, g, V_TAS, gamma_rad, rho0, B,  T0, R, V_IAS,W, h, dt)

W = W-consumed_fuel*g;
h = h + V_TAS*sin(gamma_rad)*dt;
rho = rho0*(1-B*h/T0)^(g/(R*B))*(T0/(T0-B*h));
V_TAS = V_IAS*sqrt(rho0/rho);

end



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


