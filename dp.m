
clear;
close all;
clc;

mass_1  = 0.3;                      % kg
mass_2  = 0.3;                      % kg
L1      = 50;                       % cm
L2      = 30;                       % cm
theta_1 = 30;                       % derece
theta_2 = 30;                       % derece
g       = 9.81;

% Kinematic Equations
X1 =  L1 * sin(theta_1);
Y1 = -L1 * cos(theta_1);
X2 =  L1 * sin(theta_1) + L2 * sin(theta_2);
Y2 = -L1 * cos(theta_1) - L2 * cos(theta_2);

% Velocities
X1_dot = diff(X1);
Y1_dot = diff(Y1);
X2_dot = diff(X2);
Y2_dot = diff(Y2);

% Angular Velocities
theta_1_dot = diff(theta_1);
theta_2_dot = diff(theta_2);

% Potential Energy
Energy_Pot = (-(mass_1 + mass_2) * g * L1 * cos(theta_1)) - (mass_2 * g * L2 * cos(theta_2));

% Kinetic Energy
Energy_Kin = (mass_1 * (L1^2) * (theta_1_dot^2) / 2) + ((L1^2) * (theta_1_dot^2)...
             + (L2^2) * (theta_2_dot^2) ...
             + 2 * L1 * L2 * theta_1_dot * theta_2_dot * cos(theta_1 - theta_2)) * mass_2 / 2;
         
         
 %% Lagrange Denklemleri
 
 % Parametreler
 T_Delay = 0.02;
 g     = 9.81;
 m0    = 0.5;
 m1    = 0.2;
 m2    = 0.2;
 l1    = 0.3;
 l2    = 0.3;
 L1    = 0.6;
 L2    = 0.6;
 I1    = 0.006;
 I2    = 0.006;
 t1=180;
 t2=180;
%  t0    =
%  t1    =
%  t2    =
%  t0_d  =
%  t1_d  =
%  t2_d  =
%  t0_dd =
%  t1_dd =
%  t2_dd =
 
 syms g m0 m1 m2 l1 l2 L1 L2 I1 I2 t0 t1 t2 t0_d t1_d t2_d t0_dd t1_dd t2_dd;
 
 denklem0 = ((m0+m1+m2)*t0_dd) + ((m1*l1+m2*L1)*cos(t1)*t1_dd) + m2*l2*cos(t2)*t2_dd - ((m1*l1+m2*L1)*sin(t1)*t1_d^2) - (m2*l2*sin(t2)*t2_d^2);
 denklem1 = ((m1*l1^2+m2*L1^2+I1)*t1_dd) + ((m1*l1+m2*L1)*cos(t1)*t0_dd) + ((m2*L1*l2*cos(t1-t2)*t2_dd)) + (m2*L1*l2*sin(t1-t2)*t2_d^2) - (g*(m1*l1+m2*L1)*sin(t1));
 denklem2 = (m2*l2*cos(t2)*t0_dd) + (m2*L1*l2*cos(t1-t2)*t1_dd) + ((m2*l2^2+I2)*t2_dd) - (m2*L1*l2*sin(t1-t2)*t1_d^2) - (m2*g*l2*sin(t2));
 
 Lagrange.denklem0_t0dd = solve(denklem0, t0_dd);
 Lagrange.denklem1_t1dd = solve(denklem1, t1_dd);
 Lagrange.denklem2_t2dd = solve(denklem2, t2_dd);
 
 d1 = m0 + m1 + m2;
 d2 = m1 * l1 + m2 * L1;
 d3 = m2 * l2;
 d4 = m1 * l1 * l1 + m2 * L1 * L1 + I1;
 d5 = m2 * L1 * l2;
 d6 = m2 * l2 * l2 + I2;
 f1 = (m1 * l1 + m2 * L1) * g;
 f2 = m2* l2 * g;
 
 D = [d1 d2*cos(theta1) d3*cos(theta2) ; d2*cos(theta1) d4 ...
      d5*cos(theta1-theta2) ; d3*cos(theta2) d5*cos(theta1-theta2) d6];
 C = [0 -d2*sin(theta1)*theta_1_dot -d3*sin(theta2)*theta_2_dot ; 0 0 ...
      d5*sin(theta1-theta2)*theta_2_dot ; 0 -d5*sin(theta1-theta2)*theta1_dot 0];
 G = [0 ; -f1*sin(theta1) ; -f2sin(theta2)];
 H = [1 ; 0 ; 0];
 
