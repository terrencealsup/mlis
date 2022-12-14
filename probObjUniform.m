classdef probObjUniform
    %PROBOBJUNIFORM Wrapper on uniform random variable
    %   
    
    properties
        alpha = 0; % left endpoint
        beta = 1; % right endpoint
        dim = 1;
    end
    
    methods
        function obj = probObjUniform(alpha, beta, dim)
            assert(dim == 1);
            assert(alpha < beta);
            obj.alpha = alpha;
            obj.beta = beta;
            obj.dim = dim;
        end
        
        function p = pdf(obj, x)
            assert(size(x, 2) == obj.dim)
            p = pdf('Uniform', x, obj.alpha, obj.beta);
        end
        
        function x = random(obj, m)
            x = random('Uniform', obj.alpha, obj.beta, [m, 1]);
        end
        
        function print(obj)
            disp(['alpha (left) = ', num2str(obj.alpha), ', beta (right) = ', num2str(obj.beta)]);
        end
    end
    
end

