% auxiliary function, can be ignored
clear all; close all;
addpath(genpath('FOM/'));

level = 4;
M = 1e+3;
Y = drawSamples(M);
o = zeros(M, 1);
parfor i=1:M
    disp(num2str(i));
    o(i) = beamFOM(level, Y(i, :));
end
save results/refSamples-unif Y o level M -v7.3