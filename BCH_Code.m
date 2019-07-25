m = 5;          % Parameter used to define size of GF(2^m)
t = 3;          % Number of correctable errors
k = 16;         % Number of information bits
n = 2^m - 1;    % Block length (codeword length)
p = n - k;      % Number of parity check bits
d = 2*t + 1;    % Minimum Hamming Distance between codewords

% Set up GF(2^m) that we will define the characteristics of our BCH Code %
% Primitive Polynomial = m_1(x) = x^5 + x^2 + 1 = [1 0 0 1 0 1] = 37
prim_poly = de2bi(37, 'left-msb');
M = [0 1 2 4 8 16 5 10 20 13 26 17 7 14 28 29 31 27 19 3 6 12 24 21 15 30 25 23 11 22 9 18];
GF = de2bi(M, 5, 'left-msb');

% Define the generator polynomial %
% g(x) = x^15 + x^11 + x^10 + x^9 + x^8 + x^7 + x^5 + x^3 + x^2 + x + 1 
% g(x) = 1000111110101111 (in binary) 
g_x = de2bi(36783, 'left-msb');

% ------------------------------ Encode --------------------------------- %

% Encode 16-bits of information into a 31-bit Codeword. Codeword =
% Information + Checkbits. The Checkbits are obtained by dividing a
% polynomial that represents our information by the generator polynomial,
% g(x). I will use ASCII character "A" (65, 1000001) as the information in
% this example.

% Create 16-bit binary represenation of "A"
info = de2bi(65, k, 'left-msb');

% Calculate checkbits. Append a number of zeros equal to the degree of g(x) 
% to info. This will be our dividend.
dividend = info;
dividend(end+1:end+length(g_x)-1) = 0;

% Divide by the generator polynomial, g(x) and obtain the remainder as our
% checkbits.
checkbits = polynomial_mod(dividend, g_x);

% Create final codeword by combining info & checkbits
tx_codeword = [info checkbits];

% ------------------------------ Decode --------------------------------- %

% Simulate 3 errors during the transmission of tx_codeword
rx_codeword = tx_codeword;
rx_codeword(4) = 1;
rx_codeword(9) = 1;
rx_codeword(22) = 0;

% Step 1: Compute syndrome, r'(x), from the received codeword.  For BCH, we
% are interested in components of the syndrome for error correcting. 
% Compute the 2t components of the syndrome vector, S(x), through 
% polynomial division of rx_codeword and the minimal polynomial, m_i(x), of
% each successive power of the generating element, alpha.

% minimal polynomials of alpha, alpha^2, ..., alpha^6.
m1_x = prim_poly;
m2_x = m1_x;
m4_x = m1_x;
m3_x = de2bi(61, 'left-msb');
m6_x = m3_x;
m5_x = de2bi(55, 'left-msb');

% Compute the 2t (6 in our example) components of S(x).
S1_x = polynomial_mod(rx_codeword,m1_x);
S2_x = polynomial_mod(rx_codeword,m2_x); 
S3_x = polynomial_mod(rx_codeword,m3_x);
S4_x = polynomial_mod(rx_codeword,m4_x); 
S5_x = polynomial_mod(rx_codeword,m5_x);
S6_x = polynomial_mod(rx_codeword,m6_x); 

% form system of equations in alpha (converted to Power Form)
S1_a = poly2power(GF, polynomial_mod(generateSi_a(S1_x,1), prim_poly));
S2_a = poly2power(GF, polynomial_mod(generateSi_a(S2_x,2), prim_poly));
S3_a = poly2power(GF, polynomial_mod(generateSi_a(S3_x,3), prim_poly));
S4_a = poly2power(GF, polynomial_mod(generateSi_a(S4_x,4), prim_poly));
S5_a = poly2power(GF, polynomial_mod(generateSi_a(S5_x,5), prim_poly));
S6_a = poly2power(GF, polynomial_mod(generateSi_a(S6_x,6), prim_poly));
S_a = [S1_a S2_a S3_a S4_a S5_a S6_a];


% Step 2: Find the error locator polynomial from the syndrome components 
% using Berlekamp's algorithm. 

% Initialize table.
mu = zeros(5,1);
mu(1) = -1/2;
mu(3) = 1;
mu(4) = 2;
mu(5) = 3;

sigma_x = -1.*ones(1,5,5);
sigma_x(1,5,1) = 0;
sigma_x(1,5,2) = 0;

d_mu = -1.*ones(5,1);
d_mu(1) = 0;
d_mu(2) = S1_a;

l_mu = zeros(5,1);

twomu_lmu = zeros(5,1);
twomu_lmu(1) = -1;

% Fill in the table.
for j = 2:length(mu)
	if d_mu(j) == 0
		% Calculate sigma_x(j+1).
		sigma_x(:,:,j+1) = sigma_x(:,:,j);
	else
		% Find a preceding row (row rho) with the most positive twomu_lmu
		% and d_mu ~= 0.
		rho = 0;
		most_pos = -10000;
		for k = 1:(j-1)
			if(twomu_lmu(k) > most_pos && d_mu(k) ~= -1)
				most_post = twomu_lmu(k);
				rho = k;
			end
		end

		% Calculate sigma_x(j+1).
		% Step 1: Calculate power of x^2(mu-rho)
		power = 2*(mu(j)-mu(rho));

		% Step 2: Calculate x^(power)*(sigma_x(rho))
		sigma_rho = sigma_x(:, :, rho);
		for i = 1:power
			sigma_rho(1) = sigma_rho(2);
			sigma_rho(2) = sigma_rho(3);
			sigma_rho(3) = sigma_rho(4);
			sigma_rho(4) = sigma_rho(5);
			sigma_rho(5) = -1;
		end

		% Step 3: Calculate d_mu*d_rho^(-1).
		coeff = p_mult(d_mu(j),p_inv(d_mu(rho)));
        
		% Step 4: Calculate Step3 * Step2
		for i = 1:length(sigma_rho)
			if sigma_rho(i) ~= -1
				sigma_rho(i) = p_mult(sigma_rho(i),coeff);
			end
		end
        
		% Step 5: Calculate Sigma_x(:,:,j) + Step4
		sigma_x(:,:,j+1) = poly_add(GF, sigma_x(:,:,j),sigma_rho);

	end

	% If mu == t-1, terminate the algorithm.
	if mu(j) == 2
		break; 
	end    

	% Calculate l_mu(j+1).
	degree = 0;
	for i = 1:length(sigma_x(:,:,j+1))
		if sigma_x(:,i,j+1) ~= -1
			degree = length(sigma_x(:,:,j+1)) - i;
			break;
		end
	end

	l_mu(j+1) = degree;     

	% Calculate d_mu(j+1).
	L = l_mu(j+1);
	S_sigma_coeff = -1.*ones(1,2,3);
	for i = 0:L
		S_sub = 2*mu(j) + (3-i);
		sigma_term = length(sigma_x(:,:,j+1)) - i;
		S_sigma_coeff(:,:,i+1) = [S_sub, sigma_term];
	end

	alpha_poly = -1.*ones(1,3);
	for i = 1:length(S_sigma_coeff)
		if S_sigma_coeff(1,1,i) ~= -1
			alpha_poly(i) = p_mult(S_a(S_sigma_coeff(1,1,i)), sigma_x(1,S_sigma_coeff(1,2,i),j+1));
		end                                                              
	end
    
	binary_alpha_poly = zeros(1,500);
	for i = 1:length(alpha_poly)
		if alpha_poly(i) ~= -1
			binary_alpha_poly(length(binary_alpha_poly) - alpha_poly(i)) = 1;
		end
	end

	d_mu(j+1) = poly2power(GF, polynomial_mod(binary_alpha_poly,prim_poly));

	% Calculate twomu_lmu(j+1)
	twomu_lmu(j+1) = 2*mu(j+1) - l_mu(j+1);
end

% Error locator polynomial is the last entry of sigma_x
error_loc_poly = sigma_x(:,:,length(sigma_x));


% Step 3: Use the error location polynomial to identify errant bits and
% correct them.

% Find the roots of the error locator polynomial in GF(32).
roots = -1.*ones(1,3);
count = 1;
for i = 0:length(GF)-2
	error_loc_poly_a = zeros(1, 200);
	size = length(error_loc_poly_a);
	for j = 1:length(error_loc_poly)
		if error_loc_poly(j) ~= -1
			bin_one = p_mult(error_loc_poly(j),i*(length(error_loc_poly)-j));
			error_loc_poly_a(size - bin_one) = 1;
		end
	end
	if polynomial_mod(error_loc_poly_a, prim_poly) == 0
		roots(count) = i;
		count = count+1;
	end
end

% Find inverse of roots
inv_roots = -1.*ones(1,3);
for i = 1:length(roots)
	if roots(i) ~= -1
		% Find n-Tuple Form of the root in GF
		root_i = GF(roots(i)+2, :);

		% Convert n-Tuple Form to Polynomial Form
		root_i_poly = nTuple2Poly(root_i);

		% Search through GF for root_i_poly's inverse.
		for k = 2:length(GF)
			poly = poly_mult(root_i_poly, nTuple2Poly(GF(k,:)));
			binary_poly = zeros(1,10);
			size = length(binary_poly);
			for j = 1:length(poly)
				if poly(j) ~= -1
					binary_poly(size - poly(j)) = 1;
				end
			end

			% If the mod of binary_poly and prim_poly is 1, binary_poly is
			% the multiplicative inverse.
			mod = polynomial_mod(binary_poly, prim_poly);
			inverse = 1;
			for x = 1:length(mod)
				if mod(x) == 1 && x ~= length(mod)
					inverse = 0;
				end
			end

			if inverse
				inv_roots(i) = k - 2;
				break;
			end
		end
	end
end

% Error polynomial: gives location of the errors.
e_x = inv_roots;

% XOR the received codeword and the error polynomial to fix errors!!!
binary_e_x = zeros(1,31);
for i = 1:length(e_x)
	if e_x(i) ~= -1
		binary_e_x(length(binary_e_x)- e_x(i)) = 1;
	end
end

fx_codeword = xor(rx_codeword, binary_e_x);

str1 = 'Transmitted Codeword: ';
str2 = '   Received Codeword: ';
str3 = '      Fixed Codeword: ';

disp([str1 int2str(tx_codeword)]);
disp([str2 int2str(rx_codeword)]);
disp([str3 int2str(fx_codeword)]);
