% p_mult returns the power of alpha^a*alpha^b.  
% Inputs 'a' and 'b' are powers of alpha. A value of -1 for 'a' or 'b' 
% indicates that alpha^a or alpha^b was zero in decimal, with an
% undefined power (p_mult returns -1 in this case). If a + b is a negative
% value, p_mult adds 31 until a + b is within the domain of GF(2^5).
function [ c ] = p_mult( a, b )
	if a == -1 || b == -1
		c = -1;
	elseif a == 0
		c = b;
	elseif b == 0
		c = a;
	else
		c = a + b;
		while c < 0
			c = c + 31;
	end
end
