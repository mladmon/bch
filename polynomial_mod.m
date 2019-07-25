% polynomial_mod perfoms polynomial division, a/b, and returns the modulus.
% Inputs 'a' and 'b' should be binary vectors that represent polynomials,
% where a 1 in the nth column from the right represents x^(n-1). For
% example, a = [0 1 1 0 1] represents the polynomial x^3 + x ^2 + 1. The
% function returns the remainder, 'r', of this polynomial division also in
% the form of a binary vector.
function [ r ] = polynomial_mod( a, b )
	r = a;
	while length(a) >= length(b)

		% Find the first 1 in the vector a.
		first1 = 0;
		for i = 1:length(a)
			if a(i) == 1
				first1 = i;
				break;
			end
		end

		% If 'a' contains all zeros, then we have finished the polynomial
		% division and the remainder is 0.
		if first1 == 0
			r = 0;
			break;
		end
    
		% Remove zeros past the left most one in 'a' such that the MSB of 
		% 'a' is one. NOTE: if the MSB is already 1, then 'a' remains 
		% unchanged.
		a = a(first1:length(a));

		% If removing these zeros causes 'a' to be a nonzero vector smaller
		% than 'b', then we have finished the polynomial division and the
		% remainder is 'a'.
		if length(a) < length(b)
			r = double(a);
			break;
		end

		% Pad zeros to the end of b to match dimensions of a
		extra = length(a) - length(b);
		btemp = b;
		btemp(end+1:end+extra) = 0;

		% Divide a/b by performing XOR of 'a' and 'btemp'
		a = xor(a,btemp);
	end
end
