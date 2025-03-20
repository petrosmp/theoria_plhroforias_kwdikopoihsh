clc; clear; close all;

p = 0.1;
n = 50;
samples_per_power_of_ten = 10000;
xstart = int64(0);
xend = int64(12e14);
EPSILON = 0.05;

binary_source_entropy = @(p) -p*log2(p) - (1-p)*log2(1-p);

% Start timer
start_time = tic;

cheatsheet = zeros(n+1,2);
counter = 0;
for i = 0:n
    n_choose_i = nchoosek(n, i);
    prob_i_ones = (p^i) * ((1-p)^(n-i));
    counter = counter + n_choose_i;
    cheatsheet(n-i+1, :) = [counter, prob_i_ones];
end

% Function to efficiently find probability for given x
function prob = efficient(x, cheatsheet)
    for j = 1:size(cheatsheet,1)
        if x >= cheatsheet(j,1)
            prob = cheatsheet(j,2);
            return;
        end
    end
    prob = cheatsheet(end,2);
end

x = [];
start_val = xstart;
for i = 0:log10(double(xend))
    power_of_10_increment = log10(samples_per_power_of_ten);
    until = min(10^(i+power_of_10_increment), double(xend+1));
    x = [x, start_val:int64(10^i):int64(until-1)];
    start_val = int64(until);
end

time_elapsed = toc(start_time);
disp(['Execution time: ', num2str(time_elapsed), ' s']);

figure;
ax = gca;
set(ax, 'YScale', 'log', 'XScale', 'log');
hold on;

semilogx(x, arrayfun(@(i) efficient(i, cheatsheet), x), 'b', DisplayName='probabilities');
ylabel('Probability');
xlabel('Block Index (Sorted)');
title('Probability Distribution of Block Sequences');
xlim([-10e14, 12e14]);

H = binary_source_entropy(p);
yline(2^(-n*H), 'g', DisplayName='2^{-nH}');
yline(2^(-n*(H - EPSILON)), 'r', DisplayName='2^{-n(H-\epsilon)}');
yline(2^(-n*(H + EPSILON)), 'c', DisplayName='2^{-n(H+\epsilon)}');
hold off;
legend();
