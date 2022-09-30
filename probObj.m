classdef probObj
    %PROBOBJ Probability object
    
    % example for a call
    % probCell = {{'normal', 2, {mu, sigma}}, {'lognormal', 4, {mu,sigma}}, {'gamma', 1, {alpha, beta}}
    
    properties
        probCell = {};
        indxCell = {};
        dim = 0;
    end
    
    methods
        function obj = probObj(probCell)
            obj.dim = 0;
            indx = 0;
            for i=1:length(probCell)
                name = probCell{i}{1};
                curdim = probCell{i}{2};
                param = probCell{i}{3};
                obj.indxCell{i} = indx+1:indx+curdim;
                indx = indx + curdim;
                if(strcmp(name, 'normal'))
                    probCell{i} = probObjGaussian(param{1}, param{2}, curdim);  
                elseif(strcmp(name, 'lognormal'))
                    probCell{i} = probObjLognormal(param{1}, param{2}, curdim);
                elseif(strcmp(name, 'gamma'))
                    probCell{i} = probObjGamma(param{1}, param{2}, curdim);
                elseif(strcmp(name, 'weibull'))
                    probCell{i} = probObjWeibull(param{1}, param{2}, curdim);
                elseif(strcmp(name, 'uniform'))
                    probCell{i} = probObjUniform(param{1}, param{2}, curdim);
                else
                    error('Not implemented');
                end
                obj.dim = obj.dim + curdim;
            end
            
            obj.probCell = probCell;
            
        end
        
        
        function p = pdf(obj, mu)
            assert(size(mu, 2) == obj.dim);
            
            p = ones(size(mu, 1), 1);
            for i=1:length(obj.probCell)
                p = p.*obj.probCell{i}.pdf(mu(:, obj.indxCell{i}));
            end
        end
        
        function x = random(obj, m)
            x = zeros(m, obj.dim);
            
            for i=1:length(obj.probCell)
                x(:, obj.indxCell{i}) = obj.probCell{i}.random(m);
            end
        end
        
        function print(obj)
            for i=1:length(obj.probCell)
                obj.probCell{i}.print()
            end
        end
        
        function [obj, learnStatus] = cefit(obj, x, w, I, minShrink)
            if(nargin < 5)
                minShrink = 0;
            end
            tmpProbCell = {};
            learnStatus = 1;
            for i=1:length(obj.probCell)
                [tmpProbCell{i}, curLearnStatus] = obj.probCell{i}.cefit(x(:, obj.indxCell{i}), w, I, minShrink);
                learnStatus = learnStatus & curLearnStatus;
            end
            obj.probCell = tmpProbCell;
        end
        
    end
    
end

