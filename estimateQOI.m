function qoi = estimateQOI(o, z0)
%estimateQOI estimate quantity of interest Pr(o <= z0) with samples Y.

qoi = mean(limitg(o, z0));
end

