% BME 590: Rehabilitation Robotics
% Final Project Code
% Alexander Kyu
% Date of Last Revision: 04/25/22

clear
clc
close all

%% Base Parameters

h_ua = 34; % upper arm length (cm)
m_ua = 2.9472; % upper arm weight (kg)
m_uae = 2; % upper arm exo weight (kg)
m_ua_total = m_ua + m_uae; % total mass for upper arm + exoskeleton
r1_ua = 5.5; % r1 for upper arm (cm)
r2_ua = 4.5; % r2 for upper arm (cm)
r_ua = (r1_ua + r2_ua)/2; % average r for upper arm (cm)
z_com_wo_exo_ua = 15.8704; % z center of mass for upper arm (cm)
z_com_ua = z_com_wo_exo_ua*m_ua / (m_ua_total); % z center of mass for upper arm + exoskeleton (cm)
Izz_ua = 0.5*(m_ua) * r_ua^2; % z moment of inertia (kg*cm^2)
Ixx_ua = m_ua*r_ua^2 / 4 + m_ua*h_ua^2 / 12 + m_ua*(h_ua/2)^2; % x moment of intertia (kg*cm^2)
Iyy_ua = Ixx_ua; % y moment of intertia (kg*cm^2)

h_la = 26; % lower arm length (cm)
m_la = 1.1081; % lower arm weight (kg)
m_lae = 1; % lower arm exo weight (kg)
m_la_total = m_la + m_lae; % total mass for lower arm + exoskeleton
r1_la = 4; % r1 for lower arm (cm)
r2_la = 3; % r2 for lower arm (cm)
r_la = (r1_la + r2_la)/2; % average r for lower arm (cm)
z_com_wo_exo_la = 11.7703; % z center of mass for lower arm (cm)
z_com_la = z_com_wo_exo_la*m_la / (m_la_total); % z center of mass for lower arm + exoskeleton (cm)
Izz_la = 0.5*(m_la) * r_la^2; % z moment of inertia (kg*cm^2)
Ixx_la = m_la*r_la^2 / 4 + m_la*h_la^2 / 12 + m_la*(h_la/2)^2; % x moment of intertia (kg*cm^2)
Iyy_la = Ixx_la; % y moment of intertia (kg*cm^2)

h_h = 20; % hand length (cm)
m_h = 0.55724; % hand weight (kg)
m_he = 1; % hand exo weight (kg)
m_h_total = m_h + m_he; % total mass for hand + exoskeleton
r1_h = 4.5; % r1 for hand (cm)
r2_h = 0.75; % r2 for hand (cm)
r_h = (r1_h + r2_h)/2; % average r for hand (cm)
z_com_wo_exo_h = 5.9302; % z center of mass for hand (cm)
z_com_h = z_com_wo_exo_h*m_h / (m_h_total); % z center of mass for hand + exoskeleton (cm)
Izz_h = 0.5*(m_h) * r_h^2; % z moment of inertia (kg*cm^2)
Ixx_h = m_h*r_h^2 / 4 + m_h*h_h^2 / 12 + m_h*(h_h/2)^2; % x moment of intertia (kg*cm^2)
Iyy_h = Ixx_h; % y moment of intertia (kg*cm^2)

%% Desired movement

data = xlsread('arm_spike.xlsx');
time = data(:, 1);
theta_ua = data(:, 2);
thetadot_ua = data(:, 3);
thetaddot_ua = data(:, 4);
theta_la = data(:, 5);
thetadot_la = data(:, 6);
thetaddot_la = data(:, 7);
theta_h = data(:, 8);
thetadot_h = data(:, 9);
thetaddot_h = data(:, 10);

%%  For each time point
for j = 1:length(time)

    % Parameters
    l = [h_ua/100, h_la/100, h_h/100];
    theta = [theta_ua(j), theta_la(j), theta_h(j)]; % radians
    theta_dot = [thetadot_ua(j), thetadot_la(j), thetadot_h(j)]; % radians/second
    theta_ddot = [thetaddot_ua(j), thetaddot_la(j), thetaddot_h(j)]; % radians/second^2
    d = [0, 0, 0]; % m
    d_dot = [0, 0, 0]; % m/second
    d_ddot = [0, 0, 0]; % m/s^2
    alpha = [0, 0, 0]; % radians
    a = [0, h_ua/100, h_la/100]; % m
    s = [z_com_ua/100, z_com_la/100, z_com_h/100; 0, 0, 0; 0, 0, 0]; % m 
    m = [m_ua_total, m_la_total, m_h_total]; % kg
    Ici_i(:, :, 1) = [Izz_ua/10000, 0, 0; 0, Iyy_ua/10000, 0; 0, 0, Ixx_ua/10000];
    Ici_i(:, :, 2) = [Izz_la/10000, 0, 0; 0, Iyy_la/10000, 0; 0, 0, Ixx_la/10000];
    Ici_i(:, :, 3) = [Izz_h/10000, 0, 0; 0, Iyy_h/10000, 0; 0, 0, Ixx_h/10000];
    g = 9.81; % gravity m/s^2
    sigma = [true, true, true]; % used to denote whether each joint is 
                                      % rotational(1) or prismatic(0)

    %           alpha       a       theta          d
    DH_Table = [alpha(1),   a(1),   theta(1),      d(1);
                alpha(2),   a(2),   theta(2),      d(2);
                alpha(3),   a(3),   theta(3),      d(3)];

    % Calculation of Modified DH Homogenous Transform Matrix

    for i=1:size(DH_Table,1)
        t_i = DH_Table(i, 3);
        a_i1 = DH_Table(i, 2);
        al_i1 = DH_Table(i, 1);
        d_i = DH_Table(i, 4);
        T_i1_i(:, :, i) = ...
            [cos(t_i),              -sin(t_i),             0,              a_i1;
             sin(t_i)*cos(al_i1),  cos(t_i)*cos(al_i1), -sin(al_i1),    -sin(al_i1)*d_i;
             sin(t_i)*sin(al_i1),  cos(t_i)*sin(al_i1), cos(al_i1),     cos(al_i1)*d_i;
             0,                      0,                      0,              1];
    end

    % Outward Recursion

    w0 = [0; 0; 0];
    wdot0 = [0; 0; 0];
    vdot0 = [-1*g; 0; 0];

    for i=1:size(DH_Table,1)
        R_i_i1(:, :, i) = transpose(T_i1_i(1:3, 1:3, i));
        pi1_i(:, i) = T_i1_i(1:3, 4, i);
        if i==1
            wi(:, i) = R_i_i1(:, :, i) * w0 + sigma(i)*[0; 0; theta_dot(i)];
            wdoti(:, i) = R_i_i1(:, :, i) * wdot0 + sigma(i)*R_i_i1(:, :, i) * cross(w0, [0; 0; theta_dot(i)]) + sigma(i)*[0; 0; theta_ddot(i)];
            vdoti_i(:, i) = R_i_i1(:, :, i) * (vdot0 + cross(wdot0, pi1_i(:, i)) + cross(w0, cross(w0, pi1_i(:, i)))) +... 
                                    2*~sigma(i)* R_i_i1(:, :, i)*cross(w0, [0; 0; d_dot(i)]) + ~sigma(i)*[0; 0; d_ddot(i)];
        else
            wi(:, i) = R_i_i1(:, :, i) * wi(:,i-1) + sigma(i)*[0; 0; theta_dot(i)];
            wdoti(:, i) = R_i_i1(:, :, i) * wdoti(:, i-1) + sigma(i)*R_i_i1(:, :, i) * cross(wi(:,i-1), [0; 0; theta_dot(i)]) + sigma(i)*[0; 0; theta_ddot(i)];
            vdoti_i(:, i) = R_i_i1(:, :, i) * (vdoti_i(:, i-1) + cross(wdoti(:, i-1), pi1_i(:, i)) + cross(wi(:,i-1), cross(wi(:,i-1), pi1_i(:, i)))) +...
                                   2*~sigma(i)* R_i_i1(:, :, i)*cross(wdoti(:, i-1), [0; 0; d_dot(i)]) + ~sigma(i)*[0; 0; d_ddot(i)];
        end
        vdotci_i(:, i) = vdoti_i(:, i) + cross(wdoti(:, i), s(:, i)) + cross(wi(:, i), cross(wi(:, i), s(:, i))); 
        Fi_i(:, i) = m(i)*vdotci_i(:, i);
        Ni_i(:, i) = Ici_i(:, :, i)*wdoti(:, i) + cross(wi(:, i), Ici_i(:, :, i)*wi(:, i));

    end

    % Inward Recursion

    fn1_n1 = [0; 0; 0];
    nn1_n1 = [0; 0; 0];
    pn_nt1 = [h_h/100; 0; 0];

    for i=size(DH_Table,1):-1:1
        if i == size(DH_Table,1)
            R_i_it1(:, :, i) = eye(3);
            pi_it1(:, i) = pn_nt1;
            fi_i(:, i) = R_i_it1(:, :, i)*fn1_n1 + Fi_i(:, i);
            ni_i(:, i) = R_i_it1(:, :, i)*nn1_n1 + Ni_i(:, i) + cross(pi_it1(:, i), R_i_it1(:, :, i)*fn1_n1) + cross(s(:, i), Fi_i(:, i));
        else
            R_i_it1(:, :, i) = T_i1_i(1:3, 1:3, i+1);
            pi_it1(:, i) = T_i1_i(1:3, 4, i+1);
            fi_i(:, i) = R_i_it1(:, :, i)*fi_i(:, i+1) + Fi_i(:, i);
            ni_i(:, i) = R_i_it1(:, :, i)*ni_i(:, i+1) + Ni_i(:, i) + cross(pi_it1(:, i), R_i_it1(:, :, i)*fi_i(:, i+1)) + cross(s(:, i), Fi_i(:, i));
        end
        tau(i) = sigma(i)*ni_i(3, i) + ~sigma(i)*fi_i(3, i);
    end
    
    for i=1:length(tau)
        tau_complete(j, i) = tau(i);
    end
          
end

%% Write to Files

fid = fopen('taus.txt', 'w');
for k=1:length(time)
    fprintf(fid, '%f\t%f\t%f\t%f\n', time(k), tau_complete(k, 1), tau_complete(k, 2), tau_complete(k, 3)); 
end
fclose(fid);

fid = fopen('angles.txt', 'w');
for k=1:length(time)
    fprintf(fid, '%f\t%f\t%f\t%f\n', time(k), theta_ua(k), theta_la(k), theta_h(k)); 
end
fclose(fid);

%% Plot Figures
close all

data_angles = xlsread('wm2d_angles_output.xlsx');
data_taus = xlsread('wm2d_taus_output.xlsx');
time_output = data_angles(:, 1);  % should be equivalent to time in taus_output file

angles_ua_theta = data_angles(:, 2);
angles_la_theta = data_angles(:, 3);
angles_h_theta = data_angles(:, 4);
% Didn't use thetadot from angles wm2d output because of wm2d bug of not
% correctly recording angular velocity while using angle inputs

taus_ua_theta = data_taus(:, 2);
taus_la_theta = data_taus(:, 3);
taus_h_theta = data_taus(:, 4);
taus_ua_thetadot = data_taus(:, 5);
taus_la_thetadot = data_taus(:, 6);
taus_h_thetadot = data_taus(:, 7);

figure(1)
plot(time, tau_complete(:, 1))
hold on
plot(time, tau_complete(:, 2))
plot(time, tau_complete(:, 3))
title('Joint Torques of Arm')
xlabel('Time (seconds)')
ylabel('Torque (N*m)')
legend('Upper Arm', 'Lower Arm', 'Hand', 'Location', 'best')

figure(2)
plot(time_output, angles_ua_theta)
hold on
plot(time_output, taus_ua_theta)
title('Desired vs Expected Upper Arm Joint Angle')
xlabel('Time (seconds)')
ylabel('Angle (Radians)')
legend('Desired Angle', 'Angle from Torques', 'Location', 'best')

figure(3)
plot(time_output, angles_la_theta)
hold on
plot(time_output, taus_la_theta)
title('Desired vs Expected Lower Arm Joint Angle')
xlabel('Time (seconds)')
ylabel('Angle (Radians)')
legend('Desired Angle', 'Angle from Torques', 'Location', 'best')

figure(4)
plot(time_output, angles_h_theta)
hold on
plot(time_output, taus_h_theta)
title('Desired vs Expected Wrist/Hand Joint Angle')
xlabel('Time (seconds)')
ylabel('Angle (Radians)')
legend('Desired Angle', 'Angle from Torques', 'Location', 'best')

figure(5)
plot(time, thetadot_ua)
hold on
plot(time_output, taus_ua_thetadot)
title('Desired vs Expected Upper Arm Joint Angular Velocity')
xlabel('Time (seconds)')
ylabel('Angular Velocity (Radians/Second)')
legend('Desired Angular Velocity', 'Angular Velocity from Torques', 'Location', 'best')

figure(6)
plot(time, thetadot_la)
hold on
plot(time_output, taus_la_thetadot)
title('Desired vs Expected Lower Arm Joint Angular Velocity')
xlabel('Time (seconds)')
ylabel('Angular Velocity (Radians/Second)')
legend('Desired Angular Velocity', 'Angular Velocity from Torques', 'Location', 'best')

figure(7)
plot(time, thetadot_h)
hold on
plot(time_output, taus_h_thetadot)
title('Desired vs Expected Wrist/Hand Joint Angular Velocity')
xlabel('Time (seconds)')
ylabel('Angular Velocity (Radians/Second)')
legend('Desired Angular Velocity', 'Angular Velocity from Torques', 'Location', 'best')

