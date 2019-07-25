% power2poly takes a GF(2^m) matrix and a power of alpha and returns the
% corresponding polynomial, in n-Tuple Form. GF has n-Tuple Form elements, 
% each of which are 5 bits wide.  
function [ poly ] = power2poly( GF, power )
	location = power + 2;
	poly = GF(location,:);
end
