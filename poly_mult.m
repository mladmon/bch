% poly_mult performs polynomial multiplication, a*b, correctly reduced to 
% be in the domain of GF(2)
function [ poly ] = poly_mult( a, b )
	% Expand
	poly = -1.*ones(1,25);
	for i = 1:length(a)
		for  j = 1:length(b)
			if a(i) ~= -1 && b(j) ~= -1
				poly(j+5*(i-1)) = p_mult(a(i), b(j)); 
			end
		end
	end

	% Reduce by combining like terms, subtract 2 until in domain of GF(2)
	for i = 1:length(poly)
		for  j = 1:length(poly)
			if i ~= j && poly(i) ~= -1 && poly(j) ~= -1
				if poly(i) == poly(j)
					poly(i) = -1;
					poly(j) = -1;
				end
			end
		end
	end
end
