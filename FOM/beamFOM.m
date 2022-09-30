function [ o ] = beamFOM(level, Z)
%FOM Cantilever beam FOM
% In
%   level       ...     level of discretization
%   Z = [L, t]  ...     length and thickness of beam
% Out
%   o           ...     amplitude in final oscillation

N = 2^level;
Tend = 1;
deltaT = 1e-05;
t = 0:deltaT:Tend;

%% generate system
U = ones(1, length(t));
o = zeros(size(Z, 1), 1);
for i=1:size(Z, 1)
    sys = fem_beam(N, Z(i, 1), Z(i, 2));
    y = lsim(sys, U, t);
    
    % amplitude in final oscillation
    maxima = findpeaks(y);
    if(isempty(maxima))
        maxima = max(y);
    end
    meanY = mean(y);
    o(i) = abs(maxima(end) - meanY);
end

end

