function [ poly ] = nTuple2Poly( nTuple )
%nTuple2Poly converts an n-Tuple Form entry of GF and returns the
%corresponding Polynomial Form, as a vector of powers of alpha.
   
    poly = -1.*ones(1,5);
    for j = 1:length(poly)
        if nTuple(j) == 1
            poly(j) = length(nTuple) - j;
        end
    end
end

