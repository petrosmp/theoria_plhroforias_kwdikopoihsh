syms p a
assume(p >= 0 & p <= 1);
assume(a > 0 & a < 1);

% Mutual information: I(X;Y) = H(ap) - p*H(a)
% H(x) = -x*log2(x) - (1-x)*log2(1-x)

I = -(a*p)*log(a*p)/log(2) - (1-a*p)*log(1-a*p)/log(2) + p*(a*log(a)/log(2) + (1-a)*log(1-a)/log(2));

dI_dp = diff(I, p);
p_optimal = solve(dI_dp == 0, p);

fprintf('Optimal p:\n');
pretty(p_optimal)
