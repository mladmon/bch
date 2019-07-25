% p_inv takes in a power of alpha, a, and returns the inverse of this power.
% If a = -1, then p_inv returns -1 because the power is undefined (decimal
% value was 0) and if a = 0, then p_inv returns 0 because the decimal value
% was 1. Otherwise, p_inv returns -a.
function [ b ] = p_inv( a )
	if a == -1
		b = -1;
	elseif a == 0
		b = 0;
	else
		b = -a;
end
