%{
Script with the default configuration file for the accuracy code, automatically 
called by it. It can be used as a template to overwrite some or all the default
values. All values must be stored in structure P.

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170823
%----------------------------------------------------------------------
%}

%% Directory where figures will be saved
% P.savedir = '';

%% Problem frequencies - Everything is normalized with omega = 1
P.omega_pe = 2; % electron plasma frequency
P.omega_pi = 0; % ion plasma frequency
P.omega_ce = 1.5; % electron cyclotron frequency
P.omega_ci = 0; % ion cyclotron frequency
P.nu_e = 0.1; % electron collisional frequency
P.nu_i = 0; % ion collisional frequency

%% Magnetic field direction
P.thetaB = 0*pi/180; % azimuth angle of vector B0
P.phiB = 0*pi/180; % elevation angle of vector B0

%% Unit vectors for the plane of analysis
P.u1 = [1;0;0];
P.u2 = [0;0;1];

%% Branch initial guess and parameters - Add more entries as needed
P.Bguess = cell(0);
P.n_points = [];
P.ds = [];
 
P.Bguess{end+1} = [3;0.1;P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = 0.04;
P.Bguess{end+1} = [3;0.1;P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = -0.04;

P.Bguess{end+1} = [3;0.1;pi+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = 0.04;
P.Bguess{end+1} = [3;0.1;pi+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = -0.04;

P.Bguess{end+1} = [0;1.7;pi/2+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = 0.04;
P.Bguess{end+1} = [0;1.7;pi/2+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = -0.04;

P.Bguess{end+1} = [0;1.7;-pi/2+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = 0.04;
P.Bguess{end+1} = [0;1.7;-pi/2+P.phiB];
P.n_points(end+1) = 300;
P.ds(end+1) = -0.04;

%% Numerical schemes to use
P.schemes = {};
P.h = cell(0); 
P.h{end+1} = [0.2;0.2;0.2]; 