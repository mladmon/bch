% poly2power takes a GF matrix, locates poly within this matrix, and returns
% the power corresponding to the Power Form of poly. NOTE: If poly is equal
% to [0 0 0 0 0], poly2power returns -1, which coresponds to decimal value
% of 0. GF has n-Tuple Form elements, each of which are 5 bits wide. Zeros
% must be padded onto poly, if necessary, to make poly 5 bits wide so that
% its dimensions match the elements of GF.
function [ power ] = poly2power( GF, poly )
	% Make sure poly is 5 bits wide.
	temp = bi2de(poly, 'left-msb');
	poly = de2bi(temp, 5, 'left-msb');

	% Find location of poly in GF.
	location = 0;
	for i = 1:size(GF,1)
		if GF(i,:) == poly
			location = i;
		end
	end

	power = location - 2;
end
