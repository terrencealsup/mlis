classdef probObjGaussian
    %PROBOBJGAUSSIAN Wrapper on Gaussian random variable
    %   
    
    properties
        gm = [];
        dim = 0;
    end
    
    methods
        function obj = probObjGaussian(mu, sigma, dim)
            obj.gm = gmdistribution(mu(:)', sigma, 1);
            obj.dim = dim;
        end
        
        function p = pdf(obj, x)
            assert(size(x, 2) == obj.dim)
            p = pdf(obj.gm, x);
        end
        
        function x = random(obj, m)
            x = random(obj.gm, m);
        end
        
        function print(obj)
            disp(['mu = ', num2str(obj.gm.mu), ', sigma = ', num2str(obj.gm.Sigma(:)')]);
        end
        
        function [obj, learnStatus] = cefit(obj, x, w, I, minShrink)
            if(nargin < 5)
                minShrink = 0;
            end
            
            mu = sum(bsxfun(@times, w(I == 1, :), x(I == 1, :)), 1)/sum(w(I == 1, :), 1);
            if(size(x(I == 1, :), 1) > 1)
                x = x(I == 1, :);
                w = w(I == 1, :);
            else
                % if we have a single sample only, then use all data to estimate
                % variance; better than nothing, probably
            end
            sigma = 1/sum(w)*(bsxfun(@times, w, (bsxfun(@minus, x, mu)))'*(bsxfun(@minus, x, mu)));
            % enforce symmetry in sigma to avoid numerical errors
            sigma = (sigma + sigma')/2;
            % enforce minimum sigma
            sigma = sigma + minShrink*eye(size(sigma, 1));
            disp(['      ', num2str(mu), ', ', num2str(sigma(:)')]);
            
            if(any(isnan(mu)) || any(any(isnan(sigma))) || any(any(~isreal(sigma))) || isinf(cond(sigma)))
                disp([num2str(mu), ', ', num2str(sigma(:)')]);
                gm = [];
            else
                % To cope with old Matlab versions that often fail for ill-conditioned
                % sigmas
                try
                    gm = gmdistribution(mu, sigma, size(x, 2));
                catch err
                    disp(getReport(err, 'extended'));
                    gm = [];
                end
                % try to evaluating from new gaussian model to check if Matlab thinks it is
                % ill-conditioned; if ill-conditioned use []
                try
                    pdf(gm, x(1, :));
                catch err
                    disp(getReport(err, 'extended'));
                    gm = [];
                end
                if(gm.Sigma == 0)
                    % had only one sample and therefore couldn't fix sigma
                    gm = gmdistribution(mu, obj.gm.Sigma, size(x, 2));
                end
            end
            if(isempty(gm))
                % no changes
                learnStatus = 0;
            else
                obj.gm = gm;
                learnStatus = 1;
            end
            
        end
    end
    
end

