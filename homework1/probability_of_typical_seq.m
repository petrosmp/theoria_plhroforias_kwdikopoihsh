EPSILON = 0.05;

binary_source_entropy = @(p) -p*log2(p) - (1-p)*log2(1-p); 

p = 0.45;
H = binary_source_entropy(p);

probs = [];
lower_bounds = [];
n_values = 10:10:1000;

for n = n_values
    prob = 2^-n*H;
    lower_bound = 1 - EPSILON;
    
    % Append values to vectors
    probs = [probs, prob];
    lower_bounds = [lower_bounds, lower_bound];

end

figure;
hold on;


plot(n_values, probs, 'b', DisplayName="typical seq. prob. (p=0.45)");
plot(n_values, lower_bounds, 'g', DisplayName="Lower Bound (p=0.45)");
hold off;
