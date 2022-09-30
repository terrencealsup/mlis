function [ Y, gm ] = drawSamples( M )
%DRAWSAMPLES Draw M samples from input distribution
% In
%   M           ...     number of samples
% Out
%   Y           ...     samples
%   gm          ...     distribution object


%mu = [1 0.01];
%sigma = [0.1 0; 0 0.00001];
%gm = probObj({{'normal', 2, {mu, sigma}}}); 
mu = 0.01;
sigma = 0.00001;
gm = probObj({{'weibull', 1, {1, 8}}, {'normal', 1, {mu, sigma}}, {'uniform', 1, {-2, 3}},});

Y = random(gm, M);

end

