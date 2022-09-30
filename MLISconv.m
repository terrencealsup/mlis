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



% Load the model timings
tt = load('results/timeMes.mat');

% Set to 1 to print updates, 0 otherwise
DEBUG = 1;

% Number of models
nrModels = length(modelList);

MYDEBUG('Construct models', DEBUG);

% Cell to store models
mmm = cell(length(modelList), 1);

% Check which models were selected
selI = zeros(length(modelList), 1);
for i=1:length(modelList)
    % The forward model at level modelList{i}{1}
    mmm{i} = modelList{i}{2};
    
    for j=1:length(tt.mCell)
        if(tt.mCell{j}{2}(1) == modelList{i}{1})
            selI(i, 1) = j;
        end
    end
end
MYDEBUG(['Selected models ', num2str(selI')], DEBUG);

% Get the nominal density
[~, gmNominal] = drawSamples(1);

gmCell = cell(nrModels, 1);
tCell = cell(nrModels, 1);
PfCell = cell(nrModels, 1);
PfTrueCell = cell(nrModels, 1);
VarCell = cell(nrModels, 1);
VarTrueCell = cell(nrModels, 1);
costs = zeros(nrModels, 1);
gIter = 0;

% Set the initial biasing density to the nominal density
gmPrev = gmNominal;
% Check if there is a uniform density in one of the dimensions and
% replace with Gaussian for the biasing density
for i=1:length(gmPrev.probCell)
    if isa(gmPrev.probCell{i}, 'probObjUniform')
        alpha = gmPrev.probCell{i}.alpha;
        beta = gmPrev.probCell{i}.beta;
        % Mean and variance of a uniform distribution
        mu = 0.5*(alpha + beta);
        sigma = (beta - alpha)^2/12; 
        gmPrev.probCell{i} = probObjGaussian(mu, sigma, 1);
    end
end


oldDensityFlag = 0;
tglobal = tic;

% Loop over all models
for modelIter=1:length(mmm)
    % Print time so far
    MYDEBUG(['=> ', num2str(toc(tglobal))], DEBUG);
    tglobal = tic;
    % Print current model
    MYDEBUG(['Model ', num2str(modelIter)], DEBUG);
    
    % Track iterations within each level
    stepsIter=1;
    while 1
        % Print current step within level
        MYDEBUG(['   Step ', num2str(stepsIter)], DEBUG)
        % Track global iterations
        gIter = gIter + 1;
        
        % Draw samples from current biasing density
        X = gmPrev.random(M);
        % Weight samples according to nominal density
        w = gmNominal.pdf(X)./gmPrev.pdf(X);
        % Evaluate current model at samples
        S = mmm{modelIter}(X);
        
        % Compute quantile
        if(oldDensityFlag && stepsIter > 1)
            tCell{modelIter}(stepsIter) = tCell{modelIter}(stepsIter - 2);
        else
            % Determine rho quantile of the model output samples
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
        
        % Compute indicator with current threshold
        I = limitg(S, tCell{modelIter}(stepsIter));
        
        % Compute indicator with target threshold
        Itrue = limitg(S, t);
        
        % Print number of model outputs above current threshold
        MYDEBUG(['      Have ', num2str(sum(I)), ' components left'], DEBUG);
        
        % Print number of model outputs above target threshold
        MYDEBUG(['      Have ', num2str(sum(Itrue)), ' true (w.r.t. t) components left'], DEBUG);
        
        % Get the current estimate for the threshold at the current step
        PfCell{modelIter}(stepsIter) = mean(I.*w);
        
        % Get the current estimate for the target threshold
        PfTrueCell{modelIter}(stepsIter) = mean(Itrue.*w);
        
        % Get the variance with the current threshold
        VarCell{modelIter}(stepsIter) = var(I.*w);
        
        % Get the variance with the target threshold
        VarTrueCell{modelIter}(stepsIter) = var(Itrue.*w);
        
        % Print results
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

