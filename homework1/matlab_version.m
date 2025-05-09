EPSILON = 0.05;

p = 0.1;
n = 15;

function prob = prob_of_block_with_index(i, p, n)
    binary_repr = dec2bin(i, n) - '0';
    ones_count = sum(binary_repr);
    prob = (p^ones_count) * ((1-p)^(n-ones_count));
end

function H = binary_source_entropy(p)
    H = -p*log2(p) - (1-p)*log2(1-p);
end

probabilities = [];

for i = 0:n
    probabilities = [probabilities; zeros(nchoosek(n, i), 1) + (p^i) * ((1 - p)^(n - i))];
end


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
