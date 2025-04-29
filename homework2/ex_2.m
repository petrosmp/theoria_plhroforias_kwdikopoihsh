transition_matrix = [
   1/2,  1/8,  1/8,  1/8,  1/8;
   1/4,  1/8,  1/16, 1/16, 1/2;
   1/4,  1/8,  1/8,  1/4,  1/4;
   1/8,  0,    1/2,  1/4,  1/8;
   0,    1/2,  1/4,  1/4,  0
];

stationary_distr = [1/5, 1/5, 1/5, 1/5, 1/5];

num_of_experiments = 1000;
n = 100;

function codebook = huffman_codebook(probs)
    % [probability, [symbol, code]]
    heap = cell(length(probs), 1);
    for i = 1:length(probs)
        heap{i} = {probs(i), {i, ''}};
    end
    
    [~, idx] = sort(cellfun(@(x) x{1}, heap));
    heap = heap(idx);
    
    while length(heap) > 1
        a = heap{1};
        heap(1) = [];
        b = heap{1};
        heap(1) = [];
        
        for i = 2:length(a)
            pair = a{i};
            pair{2} = ['0', pair{2}];
            a{i} = pair;
        end
        
        for i = 2:length(b)
            pair = b{i};
            pair{2} = ['1', pair{2}];
            b{i} = pair;
        end
        
        merged = {a{1} + b{1}};
        for i = 2:length(a)
            merged{end+1} = a{i};
        end
        for i = 2:length(b)
            merged{end+1} = b{i};
        end
        
        inserted = false;
        for i = 1:length(heap)
            if merged{1} <= heap{i}{1}
                heap = [heap(1:i-1); {merged}; heap(i:end)];
                inserted = true;
                break;
            end
        end
        if ~inserted
            heap{end+1} = merged;
        end
    end
    
    final_tree = heap{1};
    codes = final_tree(2:end);
    
    codebook = cell(length(probs), 1);
    for i = 1:length(codes)
        symbol = codes{i}{1};
        code = codes{i}{2};
        codebook{symbol} = code;
    end
end

function seq = generate_sequence(seqlen, transition_matrix, initial_distr)
    r = rand();
    cumul = 0;
    state = 0;
    for i = 1:length(initial_distr)
        cumul = cumul + initial_distr(i);
        if r <= cumul
            state = i;
            break;
        end
    end
    
    seq = [state];
    
    for i = 2:seqlen
        r = rand();
        cumul = 0;
        for j = 1:size(transition_matrix, 2)
            cumul = cumul + transition_matrix(state, j);
            if r <= cumul
                state = j;
                break;
            end
        end
        seq(i) = state;
    end
end

function encoded = huffman(raw, codebook)
    encoded = '';
    for i = 1:length(raw)
        encoded = [encoded, codebook{raw(i)}];
    end
end

function decoded = dehuffman(encoded, codebook)
    decoded = [];
    tmp = '';
    
    for i = 1:length(encoded)
        tmp = [tmp, encoded(i)];
        for j = 1:length(codebook)
            if strcmp(tmp, codebook{j})
                decoded = [decoded, j];
                tmp = '';
                break;
            end
        end
    end
end

function encoded = markov_huffman(raw, transition_matrix, initial_prob)
    encoded = '';
    codebook = huffman_codebook(initial_prob);
    
    for i = 1:length(raw)
        state = raw(i);
        encoded = [encoded, codebook{state}];
        codebook = huffman_codebook(transition_matrix(state, :));
    end
end

function decoded = markov_dehuffman(encoded, transition_matrix, initial_prob)
    decoded = [];
    codebook = huffman_codebook(initial_prob);
    tmp = '';
    
    i = 1;
    while i <= length(encoded)
        tmp = [tmp, encoded(i)];
        found = false;
        
        for j = 1:length(codebook)
            if strcmp(tmp, codebook{j})
                symbol = j;
                decoded = [decoded, symbol];
                tmp = '';
                codebook = huffman_codebook(transition_matrix(symbol, :));
                found = true;
                break;
            end
        end

        i = i + 1;
    end
end

stationary_codebook = huffman_codebook(stationary_distr);

stationary_avg_code_word_lengths = zeros(1, num_of_experiments);
conditional_avg_code_word_lengths = zeros(1, num_of_experiments);

stationary_entropy = 2.32;
entropy_rate = 1.875;
stationary_errors = 0;
conditional_errors = 0;

for exp = 1:num_of_experiments
    fprintf("running experiment %d\n", exp);

    initial_prob = rand(1, size(transition_matrix, 1));
    initial_prob = initial_prob / sum(initial_prob);

    seq = generate_sequence(n, transition_matrix, initial_prob);
    
    stationary_enc = huffman(seq, stationary_codebook);
    conditional_enc = markov_huffman(seq, transition_matrix, initial_prob);
    
    stationary_dec = dehuffman(stationary_enc, stationary_codebook);
    conditional_dec = markov_dehuffman(conditional_enc, transition_matrix, initial_prob);
    
    for i = 1:length(seq)
        if i <= length(stationary_dec) && seq(i) ~= stationary_dec(i)
            stationary_errors = stationary_errors + 1;
        end
        if i <= length(conditional_dec) && seq(i) ~= conditional_dec(i)
            conditional_errors = conditional_errors + 1;
        end
    end
    
    stationary_avg_code_word_lengths(exp) = length(stationary_enc) / n;
    conditional_avg_code_word_lengths(exp) = length(conditional_enc) / n;
end

fprintf('Errors using stationary huffman: %d\n', stationary_errors);
fprintf('Errors using conditional huffman: %d\n', conditional_errors);

% aggregate results for plotting
[stat_values, ~, stat_idx] = unique(stationary_avg_code_word_lengths);
stat_counts = accumarray(stat_idx, 1);
[cond_values, ~, cond_idx] = unique(conditional_avg_code_word_lengths);
cond_counts = accumarray(cond_idx, 1);

figure;
plot(stat_values, stat_counts, DisplayName='Avg Stationary Length');
hold on;
plot(cond_values, cond_counts, DisplayName='Avg Conditional Length');
legend('Avg Stationary Length', 'Avg Conditional Length');
xline(stationary_entropy, '--k', Label='Stationary Entropy', LabelOrientation='horizontal', HandleVisibility='off');
xline(entropy_rate, '--r', Label='Entropy Rate', LabelOrientation='horizontal', HandleVisibility='off');
title('Distribution of average code word length');
xlabel('Average code word length (in bits)');
ylabel('Frequency');
hold off;
