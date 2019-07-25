% poly_add takes in a GF(2^m) and two polynomials, sigma_mu and sigma_rho,
% and returns the reduced sum of these two polynomials. The two polynomials 
% must be in binary-vector form, where an 'm' in the nth column from the 
% right represents alpha^m*x^(n-1).
function [ sigma_mu_plus1 ] = poly_add( GF, sigma_mu, sigma_rho )
	sigma_mu_plus1 = -1.*ones(1,5);
	for i = 1:length(sigma_rho)
		if sigma_mu(i) == -1 && sigma_rho(i) == -1
			sigma_mu_plus1(i) = -1;
		elseif sigma_mu(i) == -1 && sigma_rho(i) ~= -1
			sigma_mu_plus1(i) = sigma_rho(i);
		elseif sigma_mu(i) ~= -1 && sigma_rho(i) == -1
			sigma_mu_plus1(i) = sigma_mu(i);
		else
			poly_mu  = power2poly(GF,sigma_mu(i));
			poly_rho = power2poly(GF,sigma_rho(i));
			poly = poly_mu + poly_rho;
			Y = poly == 2;
			poly(Y) = 0;
			sigma_mu_plus1(i) = poly2power(GF,poly);
		end
	end
end
