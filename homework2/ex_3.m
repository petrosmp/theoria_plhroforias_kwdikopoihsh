fileID = fopen('original_data_5.bin', 'rb');
A = fread(fileID, 'ubit1');
fclose(fileID);
p = 0.1;
errors = zeros(10, 2);

for n=1:10
    [hcodebook, hlens] = Huffman_codebook_generation(n, p);
    [sfcodebook, sflens] = Shannon_Fano_codebook_generation(n, p);

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
    fwrite(h_file, h_enc, 'ubit1');
    fclose(h_file);

    sf_file = fopen(sprintf('SF_compressed_data_5_%d.bin', n), 'wb');
    fwrite(sf_file, sf_enc, 'ubit1');
    fclose(sf_file);
end

fprintf("errors (one row for each n, first Huffman then S-F coding):\n");
disp(errors);

function sequences = generate_binary_sequences(n)
    total = 2^n;
    sequences = cell(1, total);
    
    for i = 0:(total-1)
        binary = dec2bin(i, n);
        sequences{i+1} = binary;
    end
end

function [code_book, code_lengths] = Shannon_Fano_codebook_generation(n, p)
    sequences = generate_binary_sequences(n);

    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    [code_book, code_lengths] = shannon_codebook(probs);
end

function [code_book, code_lengths] = Huffman_codebook_generation(n, p)
    sequences = generate_binary_sequences(n);
    
    probs = zeros(1, 2^n);
    for i = 1:2^n
        seq = sequences{i};
        ones_count = sum(seq == '1');
        probs(i) = p^ones_count * (1-p)^(n-ones_count);
    end
    
    [code_book, code_lengths] = huffman_codebook(probs);

end

function [code_book, code_lengths] = huffman_codebook(probs)
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

    code_book = zeros(length(t_codebook), max_t_codebook_len);
    code_lengths = zeros(length(t_codebook), 1);
    for i=1:length(t_codebook)
        code_word_length = strlength(strtrim(t_codebook(i, :)));
        code_book(i, 1:code_word_length) = strtrim(t_codebook(i, :))-'0';
        code_lengths(i) = code_word_length;
    end
end

function [code_book, code_lengths] = shannon_codebook(probs)
    n = length(probs);
    
    % keep the indices to restore original order
    [sorted_probs, indices] = sort(probs, 'descend');
    
    codeword_lengths = zeros(1, n);
    for i = 1:n
        codeword_lengths(i) = ceil(-log2(sorted_probs(i)));
    end
    
    codewords = cell(1, n);
    
    cum_prob = 0;
    for i = 1:n
        codeword = binary_expansion(cum_prob, codeword_lengths(i));
        codewords{i} = codeword;
        
        cum_prob = cum_prob + sorted_probs(i);
    end

    max_length = max(codeword_lengths);
    
    code_book = zeros(n, max_length);
    code_lengths = zeros(1, n);
    
    for i = 1:n
        orig_idx = indices(i);
        code_str = codewords{i};
        length_i = codeword_lengths(i);
        
        code_digits = zeros(1, max_length);
        for j = 1:length_i
            code_digits(j) = str2double(code_str(j));
        end
        
        code_book(orig_idx, :) = code_digits;
        code_lengths(orig_idx) = length_i;
    end
end

function codeword = binary_expansion(x, length)
    codeword = '';
    for i = 1:length
        x = x * 2;
        if x >= 1
            codeword = [codeword '1'];
            x = x - 1;
        else
            codeword = [codeword '0'];
        end
    end
end

function Y = encoding(X, n, code_book, code_lengths)
    % reshape the input into n-bit sequences
    xbin = reshape(X,n,[]).';
    Y = [];
    for i=1:length(xbin)
        raw = bit2int(xbin(i, :).', n);
        enc = code_book(raw+1, 1:code_lengths(raw+1));
        %fprintf("encoded %s (%d) to %s\n",  xbin(i, :)+'0', raw, enc+'0');

        if mod(i, fix(length(xbin)/100)) == 0
            fprintf("encoding sequence %d\n", i);
        end
        Y = [Y, enc];
    end
    Y = reshape(Y.', [], 1);
end

function X_decoded = decoding(Y, p, n, code_book, code_lengths)
    X_decoded = [];
    last = 1;
    i=1;
    while i <= length(Y)
        if mod(i, fix(length(Y)/100)) == 0
            fprintf("processing bit %d\n", i);
        end
        tmp = Y(last:i).';
        for j = 1:size(code_book, 1)
            if (length(tmp) == code_lengths(j)) && isequal(tmp, code_book(j, 1:code_lengths(j)))
                dec = int2bit(j-1, n);
                %fprintf("matched %s to %s, decoded to %s (%d)\n", codebook(j, 1:code_lengths(j))+'0', tmp+'0', dec+'0', j-1);
                X_decoded = [X_decoded int2bit(j-1, n)];
                last = i+1;
                break;
            end
        end
        i = i + 1;
    end
    X_decoded = reshape(X_decoded, [], 1);
end
