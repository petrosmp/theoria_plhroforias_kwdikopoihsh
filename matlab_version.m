EPSILON = 0.05;

p = 0.45;
n = 15;

% Function to compute probability of a block corresponding to index i
function prob = prob_of_block_with_index(i, p, n)
    binary_repr = dec2bin(i, n) - '0'; % Convert decimal i to binary representation
    ones_count = sum(binary_repr); % Count ones
    prob = (p^ones_count) * ((1-p)^(n-ones_count));
end

% Function to compute binary source entropy
function H = binary_source_entropy(p)
    H = -p*log2(p) - (1-p)*log2(1-p);
end

probabilities = zeros(1, 2^n);

for i = 0:(2^n - 1)
    probabilities(i+1) = prob_of_block_with_index(i, p, n);
end

probabilities = sort(probabilities, 'descend');

% Plot the results
figure;
semilogy(probabilities, 'b', DisplayName='probabilities');
xlim([-0.25e4, 3.5e4]);
ylabel('Probability');
xlabel('Block Index (Sorted)');
title('Probability Distribution of Block Sequences');

H = binary_source_entropy(p);
hold on;
yline(2^(-n*H), 'r', DisplayName='2^{-nH}');
yline(2^(-n*(H - EPSILON)), 'g', DisplayName='2^{-n(H-\epsilon)}');
yline(2^(-n*(H + EPSILON)), 'm', DisplayName='2^{-n(H+\epsilon)}');
hold off;
legend();
