fileID = fopen('original_data_5.bin', 'rb');
A = fread(fileID, 'ubit1');
fclose(fileID);
p = 0.1;
errors = zeros(10, 2);

for n=1:10
    [hcodebook, hlens] = huffman_codebook_seq(n, p);
    [sfcodebook, sflens] = shannon_fano_codebook_seq(n, p);

    h_enc = encoding(A, n, hcodebook, hlens);
    % h_dec = decoding(h_enc, p, n, hcodebook, hlens);
    % 
    % for i = 1:length(A)
    %     if h_dec(i) ~= A(i)
    %         errors(n, 1) = errors(n, 1) + 1;
    %     end
    % end

    sf_enc = encoding(A, n, sfcodebook, sflens);
    % sf_dec = decoding(sf_enc, p, n, sfcodebook, sflens);
    % 
    % for i = 1:length(A)
    %     if sf_dec(i) ~= A(i)
    %         errors(n, 2) = errors(n, 2) + 1;
    %     end
    % end

    h_file = fopen(sprintf('H_compressed_data_5_%d.bin', n), 'wb');
    fwrite(h_file, h_enc, 'ubit1');
    fclose(h_file);

    sf_file = fopen(sprintf('SF_compressed_data_5_%d.bin', n), 'wb');
    fwrite(sf_file, sf_enc, 'ubit1');
    fclose(sf_file);

end

fprintf("errors (one row for each n, first Huffman then S-F coding):\n");
disp(errors);

function [left, right] = split_group(group)
    % split a group into two subgroups with nearly equal probability
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
    
    % add 0 to left group
    for i = 1:length(left)
        left{i}{3} = [left{i}{3}, 0];
    end
    
    % add 1 to right group
    for i = 1:length(right)
        right{i}{3} = [right{i}{3}, 1];
    end
end

function [codebook, lengths] = shannon_fano_codebook(probs)
    [sorted_probs, indices] = sort(probs, 'descend');
    
    tuples = cell(1, length(probs));
    for i = 1:length(probs)
        tuples{i} = {indices(i), sorted_probs(i), ''};
    end
    
    groups = {tuples};
    t_codebook = cell(1, length(probs));
    
    while ~isempty(groups)
        group = groups{end};
        groups(end) = [];
        
        [left, right] = split_group(group);
        
        if isscalar(left)
            t_codebook{left{1}{1}} = left{1}{3};
        else
            groups{end+1} = left;
        end
        
        if isscalar(right)
            t_codebook{right{1}{1}} = right{1}{3};
        else
            groups{end+1} = right;
        end
    end
    
    t_codebook = string(t_codebook).';

    max_t_codebook_len = max(cellfun(@(x) strlength(x), t_codebook));

    t_codebook = char(cellstr(t_codebook));

    codebook = zeros(length(t_codebook), max_t_codebook_len);
    lengths = zeros(length(t_codebook), 1);
    for i=1:length(t_codebook)
        code_word_length = strlength(strtrim(t_codebook(i, :)));
        codebook(i, 1:code_word_length) = strtrim(t_codebook(i, :));
        lengths(i) = code_word_length;
    end
end

function [codebook, lengths] = shannon_fano_codebook_seq(n, p)
    sequences = generate_binary_sequences(n);

    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    [codebook, lengths] = shannon_fano_codebook(probs);
end

function sequences = generate_binary_sequences(n)
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
    
    t_codebook = cell(length(probs), 1);
    for i = 1:length(codes)
        symbol = codes{i}{1};
        code = codes{i}{2};
        t_codebook{symbol} = code;
    end


    
    t_codebook = string(t_codebook).';

    max_t_codebook_len = max(cellfun(@(x) strlength(x), t_codebook));

    t_codebook = char(cellstr(t_codebook));

    codebook = zeros(length(t_codebook), max_t_codebook_len);
    lengths = zeros(length(t_codebook), 1);
    for i=1:length(t_codebook)
        code_word_length = strlength(strtrim(t_codebook(i, :)));
        codebook(i, 1:code_word_length) = strtrim(t_codebook(i, :))-'0';
        lengths(i) = code_word_length;
    end
end

function [codebook, lengths] = huffman_codebook_seq(n, p)
    sequences = generate_binary_sequences(n);
    
    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    [codebook, lengths] = huffman_codebook(probs);

end

function code = encoding(X, n, codebook, code_lengths)
    % reshape the input into n-bit sequences
    xbin = reshape(X,n,[]).';
    code = [];
    for i=1:length(xbin)
        raw = bit2int(xbin(i, :).', n);
        enc = codebook(raw+1, 1:code_lengths(raw+1));
        %fprintf("encoded %s (%d) to %s\n",  xbin(i, :)+'0', raw, enc+'0');

        if mod(i, fix(length(xbin)/100)) == 0
            fprintf("encoding sequence %d\n", i);
        end
        code = [code, enc];
    end
    code = reshape(code.', [], 1);
end

function decoded = decoding(Y, p, n, codebook, code_lengths)
    decoded = [];
    last = 1;
    i=1;
    while i <= length(Y)
        if mod(i, fix(length(Y)/100)) == 0
            fprintf("processing bit %d\n", i);
        end
        tmp = Y(last:i).';
        for j = 1:size(codebook, 1)
            if (length(tmp) == code_lengths(j)) && isequal(tmp, codebook(j, 1:code_lengths(j)))
                dec = int2bit(j-1, n);
                %fprintf("matched %s to %s, decoded to %s (%d)\n", codebook(j, 1:code_lengths(j))+'0', tmp+'0', dec+'0', j-1);
                decoded = [decoded int2bit(j-1, n)];
                last = i+1;
                break;
            end
        end
        i = i + 1;
    end
    decoded = reshape(decoded, [], 1);
end
