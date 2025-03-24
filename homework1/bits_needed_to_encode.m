EPSILON = 0.05;

% Binary source entropy function
binary_source_entropy = @(p) -p*log2(p) - (1-p)*log2(1-p); 

% Set probability value
p = 0.45;
H = binary_source_entropy(p);

% Initialize vectors to store probabilities and lower bounds
bits_needed = [];
lower_bounds = [];
n_values = 10:10:1000;

% Loop from n = 10 to n = 1000 with a step size of 10
for n = n_values
    bits = n*H;  % Probability calculation
    lower_bound = 1 - EPSILON;  % Lower bound
    
    % Append values to vectors
    bits_needed = [bits_needed, bits];

end

figure;
hold on;


plot(n_values, bits_needed, 'b', DisplayName="typical seq. prob. (p=0.45)");
yline(H, 'c', DisplayName='H');
yline(H - EPSILON, 'c', DisplayName='H-\epsilon');
yline(H + EPSILON, 'c', DisplayName='H+\epsilon');
hold off;
