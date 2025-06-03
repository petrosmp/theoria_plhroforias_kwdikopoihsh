rates = [1, 1/2, 1/3, 1/4, 1/5, 1/8];
p = 1/10;
num_of_experiments = 100000;
bit_error_rates = zeros(length(rates), 10);

for rate_idx = 1:length(rates)
    rate = rates(rate_idx);
    unit = floor(1 / rate);

    for n = unit:unit:10*unit
        k = floor(n*rate);
        disp(['rate ', num2str(rate), ', n ', num2str(n), ', unit ', num2str(unit), ...
              ', (', num2str(2^k), ') codewords']);

        n_r_bit_errors = 0;
        parfor i = 1:num_of_experiments            % can be replaced with parfor for faster execution
            random_code = random_N_R_code(n, rate);
            experiment_bit_errors = 0;
            for j = 1:10

                message = randi([0 1], k, 1);                           % generate random k-bit message
                cw_idx = bit2int(message, k)+1;                         % make it decimal so it can be used as an index
                cw = random_code(cw_idx, :);                            % get the codeword of that index (i.e. of that message)
                channeled = channel(cw, p);                             % transmit it
                msg_estimate = ML_decode(channeled, random_code, k);    % estimate the message based on the noisy codeword with ML 

                bit_errors = sum(message ~= msg_estimate);
                experiment_bit_errors = experiment_bit_errors + bit_errors;
            end
            n_r_bit_errors = n_r_bit_errors + experiment_bit_errors;
        end
        bit_error_rates(rate_idx, n/unit) = bit_error_rates(rate_idx, n/unit) + n_r_bit_errors/(n*num_of_experiments*10);
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

function estimate = ML_decode(received, inputs, k)
    distances = sum(received ~= inputs, 2);
    [~, min_idx] = min(distances);
    estimate = int2bit(min_idx+1, k); % make the estimated codeword index back
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

