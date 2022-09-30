classdef probObjWeibull
    %PROBOBJGAUSSIAN Wrapper on Gaussian random variable
    %   
    
    properties
        alpha = 0; % scale
        beta = 0; % shape
        dim = 0;
    end
    
    methods
        function obj = probObjWeibull(alpha, beta, dim)
            assert(dim == 1);
            obj.alpha = alpha;
            obj.beta = beta;
            obj.dim = dim;
        end
        
        function p = pdf(obj, x)
            assert(size(x, 2) == obj.dim)
            p = pdf('Weibull', x, obj.alpha, obj.beta);
        end
        
        function x = random(obj, m)
            x = random('Weibull', obj.alpha, obj.beta, [m, 1]);
        end
        
        function print(obj)
            disp(['alpha (scale) = ', num2str(obj.alpha), ', beta (shape) = ', num2str(obj.beta)]);
        end
        
        function [obj, learnStatus] = cefit(obj, x, w, I, minShrink)
            parmHat = wblfit(x(I == 1, :), [], [], w(I == 1, :));
            obj.alpha = parmHat(1);
            obj.beta = parmHat(2);
            %obj.beta = min(parmHat(2), 10);
            learnStatus = 1;
        end
    end
    
end

