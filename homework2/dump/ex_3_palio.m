fileID = fopen('original_data_5.bin');
A = fread(fileID, 'ubit1');
p = 0.1;
errors = zeros(10, 2);


for n=1:10
    [hcodebook, hlens] = huffman_codebook_seq(n, p);
    [sfcodebook, sflens] = shannon_fano_codebook_seq(n, p);

    h_enc = encoding(A, n, hcodebook, hlens);
    h_dec = decoding(h_enc, p, n, hcodebook, hlens);

    for i = 1:length(A)
        if h_dec(i) ~= A(i)
            errors(n, 1) = errors(n, 1) + 1;
        end
    end

    sf_enc = encoding(A, n, sfcodebook, sflens);
    sf_dec = decoding(sf_enc, p, n, sfcodebook, sflens);

    for i = 1:length(A)
        if sf_dec(i) ~= A(i)
            errors(n, 2) = errors(n, 2) + 1;
        end
    end

    h_file = fopen(sprintf('H_compressed_data_5_%d.bin', n), 'wb');
    fwrite(h_file, h_enc-'0', 'ubit1');
    fclose(h_file);

    sf_file = fopen(sprintf('SF_compressed_data_5_%d.bin', n), 'wb');
    fwrite(sf_file, sf_enc-'0', 'ubit1');
    fclose(sf_file);

end

function entropy = binary_source_entropy(p)
    % Calculate binary entropy function
    entropy = -p * log2(p) - (1-p) * log2(1-p);
end

function [left, right] = split_group(group)
    % Split a group into two subgroups with nearly equal probability
    best_diff = inf;
    split_index = 0;
    total = sum(cellfun(@(x) x{2}, group));
    running_sum = 0;
    
    for i = 1:(length(group)-1)
        running_sum = running_sum + group{i}{2};
        diff = abs(running_sum - (total - running_sum));
        
        if diff < best_diff
            best_diff = diff;
            split_index = i;
        end
    end
    
    left = group(1:split_index);
    right = group(split_index+1:end);
    
    % Add 0 to left codes
    for i = 1:length(left)
        left{i}{3} = [left{i}{3}, '0'];
    end
    
    % Add 1 to right codes
    for i = 1:length(right)
        right{i}{3} = [right{i}{3}, '1'];
    end
end

function [codebook, lengths] = shannon_fano_codebook(probs)
    % Create Shannon-Fano codebook from probability distribution
    % Sort probabilities in descending order
    [sorted_probs, indices] = sort(probs, 'descend');
    
    % Initialize tuples
    tuples = cell(1, length(probs));
    for i = 1:length(probs)
        tuples{i} = {indices(i), sorted_probs(i), ''};
    end
    
    % Initialize groups and codebook
    groups = {tuples};
    codebook = cell(1, length(probs));
    
    while ~isempty(groups)
        group = groups{end};
        groups(end) = [];
        
        [left, right] = split_group(group);
        
        if length(left) == 1
            codebook{left{1}{1}} = left{1}{3};
        else
            groups{end+1} = left;
        end
        
        if length(right) == 1
            codebook{right{1}{1}} = right{1}{3};
        else
            groups{end+1} = right;
        end
    end
    
    codebook = string(codebook).';
    codebook = char(cellstr(codebook));
    lengths = zeros(length(codebook), 1);
    for i=1:length(codebook)
        lengths(i) = strlength(strtrim(codebook(i, :)));
    end

end

function [codebook, lengths] = shannon_fano_codebook_seq(n, p)
    % Generate Shannon-Fano codebook for binary sequences of length n
    % with probability p of 1
    
    % Generate all possible binary sequences
    sequences = generate_binary_sequences(n);
    
    % Calculate probabilities
    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    % Create the codebook
    [codebook, lengths] = shannon_fano_codebook(probs);
end

function sequences = generate_binary_sequences(n)
    % Generate all binary sequences of length n
    total = 2^n;
    sequences = cell(1, total);
    
    for i = 0:(total-1)
        binary = dec2bin(i, n);
        sequences{i+1} = binary;
    end
end

function [codebook, lengths] = huffman_codebook(probs)
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

    codebook = string(codebook);
    codebook = char(cellstr(codebook));

    lengths = zeros(length(codebook), 1);
    for i=1:length(codebook)
        lengths(i) = strlength(strtrim(codebook(i, :)));
    end

end

function [codebook, lengths] = huffman_codebook_seq(n, p)
    % Generate Huffman codebook for binary sequences of length n
    % with probability p of 1
    
    % Generate all possible binary sequences
    sequences = generate_binary_sequences(n);
    
    % Calculate probabilities
    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    % Create the codebook
    [codebook, lengths] = huffman_codebook(probs);

end

function code = encoding(X, n, codebook, code_lengths)
    % reshape the input into n-bit sequences
    xbin = reshape(X,n,[]).';
    code = [];
    for i=1:length(xbin)
        raw = bit2int(xbin(i, :).', n);
        enc = codebook(raw+1, :);
        %fprintf("encoded %s (%d) to %s\n",  char(xbin(i, :)+'0'), raw, codebook(raw+1, :));

        if mod(i, fix(length(xbin)/100)) == 0
            fprintf("encoding sequence %d\n", i);
        end
        code = [code, strtrim(enc)];
    end
    code = reshape(code.', [], 1);
end

function decoded = decoding(Y, p, n, codebook, code_lengths)
    % Decode a codeword using the given codebook
    decoded = [];
    last = 1;
    i=1;
    while i <= length(Y)
        if mod(i, fix(length(Y)/100)) == 0
            fprintf("processing bit %d\n", i);
        end
        tmp = Y(last:i).';
        for j = 1:length(codebook)
            if strcmp(strtrim(codebook(j, :)), tmp) % this is extremely slow, do it with non chars
                %fprintf("matched %s to %s\n", codebook(j), tmp);
                %fprintf("decoded %s to %s %d\n", tmp, char(int2bit(j-1, n))+'0', j-1)
                decoded = [decoded int2bit(j-1, n)];
                last = i+1;
                break;
            end
        end
        i = i + 1;
    end
    decoded = reshape(decoded, [], 1);
end
