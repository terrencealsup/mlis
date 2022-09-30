classdef probObjLognormal
    %PROBOBJGAUSSIAN Wrapper on Gaussian random variable
    %   
    
    properties
        gm = [];
    end
    
    methods
        function obj = probObjLognormal(mu, sigma, dim)
           obj.gm = probObjGaussian(mu, sigma, dim);
        end
        
        function p = pdf(obj, x)
            assert(size(x, 2) == obj.gm.dim)
            p = obj.gm.pdf(log(x))./prod(x, 2);
        end
        
        function x = random(obj, m)
            x = random(obj.gm, m);
            x = exp(x);
        end
        
        function print(obj)
            obj.gm.print()
        end
        
        function [obj, learnStatus] = cefit(obj, x, w, I, minShrink)
            if(nargin < 5)
                minShrink = 0;
            end
            [obj.gm, learnStatus] = obj.gm.cefit(log(x), w, I, minShrink);
        end
    end
    
end

