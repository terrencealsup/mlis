function [ yesno ] = limitg( o, z0 )
%LIMITG Checks if o <= z0


yesno = zeros(size(o));
yesno(o <= z0) = 1;

end

