rates = [1, 1/2, 1/3, 1/4, 1/5, 1/8];
p = 1/10;
num_of_experiments = 100000;
bit_error_rates = zeros(length(rates), 10);

for rate_idx = 1:length(rates)
    rate = rates(rate_idx);
    unit = floor(1 / rate);

    for n = unit:unit:10*unit
        codewords = random_N_R_code_distinct(n, rate);
        num_codewords = 2^floor(n*rate);
        disp(['rate ', num2str(rate), ', n ', num2str(n), ', unit ', num2str(unit), ...
              ', ', num2str(size(codewords,1)), ' (', num2str(num_codewords), ...
              ') codewords, ', num2str(length(unique(codewords, 'rows'))), ' distinct']);

        for i = 1:num_of_experiments
            for j = 1:10
                cw_idx = randi(size(codewords,1));
                cw = codewords(cw_idx, :);
                channeled = channel(cw, p);
                decoded = ML_decode(channeled, codewords);
                bit_errors = sum(cw ~= decoded);
                bit_error_rates(rate_idx, n/unit) = bit_error_rates(rate_idx, n/unit) + bit_errors/(n*num_of_experiments*10);
            end
        end
    end
    disp(['bit error rates for rate ', num2str(rate), ': ',  mat2str(bit_error_rates(rate_idx, :))]);
end


function codewords = random_N_R_code_distinct(n, r)
    num_codewords = 2^floor(n*r);
    codeword_set = zeros(num_codewords, n);
    picked = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    i = 1;
    while i <= num_codewords
        bits = randi([0, 1], 1, n);
        bin_str = char(bits + '0');
        key = string(bin_str);
        if ~isKey(picked, key)
            picked(key) = true;
            codeword_set(i, :) = bits;
            i = i + 1;
        end
    end
    codewords = codeword_set;
end

function codewords = random_N_R_code(n, r)
    num_codewords = 2^floor(n*r);
    codewords = randi([0,1], num_codewords, n);
end

function output = channel(x, p)
    flip = rand(1, length(x)) < p;
    output = xor(x, flip);
end

function estimate = ML_decode(received, inputs)
    min_dist = sum(received ~= inputs(1,:));
    estimate = inputs(1,:);
    for i = 2:size(inputs,1)
        dist = sum(received ~= inputs(i,:));
        if dist < min_dist
            min_dist = dist;
            estimate = inputs(i,:);
        end
    end
end


figure;
hold on;

for i = 1:length(rates)
    plot(1:10, bit_error_rates(i, :), DisplayName=['R = ' num2str(rates(i), 3)], LineWidth=1.5);
end

title('Evolution of Bit Error Rate with Codeword Length');
xlabel('Sequence Length (in bits, xR^{-1})');
ylabel('Bit Error Rate');
legend('show');
grid on;
hold off;

