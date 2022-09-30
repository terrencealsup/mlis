classdef probObjGamma
    %PROBOBJGAUSSIAN Wrapper on Gaussian random variable
    %   
    
    properties
        alpha = 0;
        beta = 0;
        dim = 0;
    end
    
    methods
        function obj = probObjGamma(alpha, beta, dim)
            assert(dim == 1);
            obj.alpha = alpha;
            obj.beta = beta;
            obj.dim = dim;
        end
        
        function p = pdf(obj, x)
            assert(size(x, 2) == obj.dim)
            p = pdf('Gamma', x, obj.alpha, obj.beta);
        end
        
        function x = random(obj, m)
            x = random('Gamma', obj.alpha, obj.beta, [m, 1]);
        end
        
        function print(obj)
            disp(['alpha = ', num2str(obj.alpha), ', beta = ', num2str(obj.beta)]);
        end
        
        function [obj, learnStatus] = cefit(obj, x, w, I, minShrink)
            % we only search for beta, the alpha parameter is fixed
            obj.beta = sum(bsxfun(@times, w(I == 1, :), x(I == 1, :)), 1)/(sum(w(I == 1, :), 1)*obj.alpha);
            learnStatus = 1;
        end
    end
    
end

