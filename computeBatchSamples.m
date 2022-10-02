function [] = computeBatchSamples(jobnum)
% computeBatchSamples draws a batch of samples

addpath(genpath('FOM/'));

rng(jobnum);

level = 4;
M = 1e+3;
Y = drawSamples(M);
o = zeros(M, 1);
for i=1:M
    disp(num2str(i));
    o(i) = beamFOM(level, Y(i, :));
end
file_name = sprintf('results/batchSamples/batchSamples-%d.mat', jobnum);
save(file_name, 'Y', 'o', 'level', 'M', '-v7.3');

end

