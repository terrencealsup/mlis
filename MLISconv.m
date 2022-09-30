function [ PfCell, VarCell, gmCell, costs, PfTrueCell, VarTrueCell, tCell ] = MLISconv( modelList, M, t, delta, rho, sigmaFactor, minShrink )
%TESTMLIS Multifidelity Importance Sampling
% In
%   modelList       ...     list of models
%   M               ...     number of samples on each level
%   t               ...     threshold parameter for failure prob
%   delta           ...     delta parameter (min progress)
%   rho             ...     (1-rho) quantiles used as step
%   sigmaFactor     ...     increase covariance of biasing by sigmaFactor
%   minShrink       ...     minimal sigma
% Out
%   PfCell          ...     list of failure probability on levels
%   VarCell         ...     list of variances on levels
%   gmCell          ...     cell with biasing densities
%   costs           ...     costs on each level
%   PfTrueCell      ...     list of est fail prob with t
%   VarTrueCell     ....    list of var with t

%addpath(genpath('~/uni/apro/mlib'));
%addpath(genpath('FOM/'));
%addpath('../code');
tt = load('results/timeMes.mat');

%% parameters and models
DEBUG = 1;
nrModels = length(modelList);
% create models
MYDEBUG('Construct models', DEBUG);
mmm = cell(length(modelList), 1); %constructModels(modelList);
%selI = selectModels(tt.mCell, modelList);
selI = zeros(length(modelList), 1);
for i=1:length(modelList)
    mmm{i} = modelList{i}{2};
    for j=1:length(tt.mCell)
        if(tt.mCell{j}{2}(1) == modelList{i}{1})
            selI(i, 1) = j;
        end
    end
end
MYDEBUG(['Selected models ', num2str(selI')], DEBUG);

%% construct steps in CE
[~, gmNominal] = drawSamples(1);
gmCell = cell(nrModels, 1);
tCell = cell(nrModels, 1);
PfCell = cell(nrModels, 1);
PfTrueCell = cell(nrModels, 1);
VarCell = cell(nrModels, 1);
VarTrueCell = cell(nrModels, 1);
costs = zeros(nrModels, 1);
gIter = 0;
gmPrev = gmNominal;
oldDensityFlag = 0;
tglobal = tic;
for modelIter=1:length(mmm)
    MYDEBUG(['=> ', num2str(toc(tglobal))], DEBUG);
    tglobal = tic;
    MYDEBUG(['Model ', num2str(modelIter)], DEBUG);
    stepsIter=1;
    while 1
        MYDEBUG(['   Step ', num2str(stepsIter)], DEBUG)
        gIter = gIter + 1; % global iterations
        
        X = gmPrev.random(M); % draw samples
        w = gmNominal.pdf(X)./gmPrev.pdf(X); % weights
        S = mmm{modelIter}(X); % evaluate model
        
        % compute quantile
        if(oldDensityFlag && stepsIter > 1)
            tCell{modelIter}(stepsIter) = tCell{modelIter}(stepsIter - 2);
        else
            tCell{modelIter}(stepsIter) = quantile(S, rho);
        end
        
        if(stepsIter > 1 && ~oldDensityFlag) % check if we made progress
            minTsize = tCell{modelIter}(stepsIter - 1) - delta;
            if(isnan(tCell{modelIter}(stepsIter)))
                MYDEBUG(['      Proposed ', num2str(tCell{modelIter}(stepsIter)), ', reset to ', num2str(minTsize)], DEBUG);
                tCell{modelIter}(stepsIter) = min(tCell{modelIter}(stepsIter), minTsize);
                error('new t is nan');
            end
            if(tCell{modelIter}(stepsIter) > minTsize)
                MYDEBUG(['      Proposed ', num2str(tCell{modelIter}(stepsIter)), ' but larger than ', num2str(minTsize)], DEBUG);
                %MYDEBUG(['      Abort'], DEBUG);
                %error('Would require stepping forward with delta, better to abort');
                tCell{modelIter}(stepsIter) = min(tCell{modelIter}(stepsIter), minTsize);
            end
        end
        if(t > tCell{modelIter}(stepsIter))
            MYDEBUG(['      Proposed ', num2str(tCell{modelIter}(stepsIter)), ' but below ', num2str(t)], DEBUG);
            
            OneMinusRhoQuantile = quantile(S, 1-rho);
            if(t < OneMinusRhoQuantile)
                MYDEBUG(['      Between rho, ', num2str(tCell{modelIter}(stepsIter)), ', and 1-rho quantile, ', num2str(OneMinusRhoQuantile)], DEBUG);
                tCell{modelIter}(stepsIter) = t;
            else
                tCell{modelIter}(stepsIter) = OneMinusRhoQuantile;
                MYDEBUG(['      Using 1-rho quantile ', num2str(tCell{modelIter}(stepsIter))], DEBUG);
            end
        end
        MYDEBUG(['      Selected tCell{', num2str(modelIter), '}(', num2str(stepsIter), ') = ', num2str(tCell{modelIter}(stepsIter))], DEBUG);
        % reset oldDensityFlag
        oldDensityFlag = 0;
        
        % compute indicator
        I = limitg(S, tCell{modelIter}(stepsIter));
        Itrue = limitg(S, t);
        MYDEBUG(['      Have ', num2str(sum(I)), ' components left'], DEBUG);
        MYDEBUG(['      Have ', num2str(sum(Itrue)), ' true (w.r.t. t) components left'], DEBUG);
        
        % derive estimate on this level
        PfCell{modelIter}(stepsIter) = mean(I.*w);
        PfTrueCell{modelIter}(stepsIter) = mean(Itrue.*w);
        VarCell{modelIter}(stepsIter) = var(I.*w);
        VarTrueCell{modelIter}(stepsIter) = var(Itrue.*w);
        MYDEBUG(['      Pf = ', num2str(PfCell{modelIter}(stepsIter)), ...
            ', Var = ', num2str(VarCell{modelIter}(stepsIter))], DEBUG);
        MYDEBUG(['      PfTrue = ', num2str(PfTrueCell{modelIter}(stepsIter)), ...
            ', VarTrue = ', num2str(VarTrueCell{modelIter}(stepsIter))], DEBUG);
        MYDEBUG(['      est quantile = ', num2str(mean(I)), ...
            ', rho = ', num2str(rho)], DEBUG);
        
        % learn new density
        [gmCell{modelIter, stepsIter}, learnStatus] = gmPrev.cefit(X, w, I, minShrink);
        if(learnStatus == 0)
            disp('      No density learned, reverting to previous density');
            gmCell{modelIter, stepsIter} = gmPrev;
            oldDensityFlag = 1; % set flag that we do NOT enforce delta in next step
%         elseif(gmCell{modelIter, stepsIter}.Sigma == 0)
%             % had only one sample and therefore couldn't fix sigma
%             gmCell{modelIter, stepsIter} = gmdistribution(gmCell{modelIter, stepsIter}.mu, gmPrev.Sigma, size(X, 2));
        end
        if(sigmaFactor ~= 1)
            error('Sigma factor is not supported anymore');
        end
%         % increase sigma by factor sigmaFactor, to increase robustness
%         MYDEBUG(['      Apply sigma factor ', num2str(sigmaFactor)], DEBUG);
%         gmCell{modelIter, stepsIter} = gmdistribution(gmCell{modelIter, stepsIter}.mu, sigmaFactor*gmCell{modelIter, stepsIter}.Sigma, size(X, 2));
%         %gmCell{modelIter, stepsIter} = gmdistribution(gmCell{modelIter, stepsIter}.mu, gmPrev.Sigma, size(X, 2));
        % plot learned mu and sigma
        if(DEBUG)
            gmCell{modelIter, stepsIter}.print()
        end
        %MYDEBUG(['      Learned density mu ', num2str(gmCell{modelIter, stepsIter}.mu), ', sigma ', num2str(gmCell{modelIter, stepsIter}.Sigma(:)')], DEBUG);
        
        
        % bookkeeping
        costs(gIter) = tt.tt(selI(modelIter))*M;
        gmPrev = gmCell{modelIter, stepsIter};
             
        
        % check if we reached t
        if(tCell{modelIter}(stepsIter) == t)
            break
        end
        stepsIter = stepsIter + 1;
    end
end


end

