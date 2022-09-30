% auxiliary function for measuring costs of each model; ignore
clear variables; close all;

addpath(genpath('FOM/'));

runs = 10;
levelList = [2, 3, 4];

yMat = zeros(runs, length(levelList));
mCell = cell(length(levelList), 1);
tt = zeros(length(levelList), 1);
tterr = zeros(length(levelList), 1);
mIndx = 1;

Y = drawSamples(runs);

%% Coarse grid approximations
for levelIter=1:length(levelList)
    level = levelList(levelIter);
    disp(num2str(level));
    
    % compute samples and measure time
    t = tic;
    yMat(:, mIndx) = beamFOM(level, Y);
    tt(mIndx) = toc(t)/runs;
    tterr(mIndx) = norm(yMat(:, 1) - yMat(:, mIndx), 2)/norm(yMat(:, 1), 2);
    mCell{mIndx} = {'CoarseGrid', [level, 0]};
    mIndx = mIndx + 1;
end

save results/timeMes.mat yMat tt tterr mCell levelList;
