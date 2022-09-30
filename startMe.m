clear all; close all;
addpath(genpath('FOM/'));

modelList = {{2, @(z)beamFOM(2, z)}, {3, @(z)beamFOM(3, z)}, {4, @(z)beamFOM(4, z)}};

M = 500; % number of samples per level
t = 4e-07; % threshold, i.e., find P[f(x) < t]

% the following are most likely not tuning parameters
rho = 0.1; % (1-rho) quantiles used as step which means that M has to be large enough to faithfully estimate an event with probability rho; if approach fails, slightly increase rho to, e.g., 0.2, take more steps per level
delta = 1e-02; % minimal step size
sigmaFactor = 1; % overestimate the sigma for the biasing density by this much 
minShrink = 1e-3; % the minimal covariance the biasing density can have; prevents collapsing to a point 

% rough reference: 4.4000e-04 for t = 4e-07
PfCell = MLISconv(modelList, M, t, delta, rho, sigmaFactor, minShrink);
