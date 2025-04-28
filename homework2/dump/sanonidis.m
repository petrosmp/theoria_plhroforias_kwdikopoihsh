
function sequences = generate_binary_sequences(n)
    total = 2^n;
    sequences = cell(1, total);
    
    for i = 0:(total-1)
        binary = dec2bin(i, n);
        sequences{i+1} = binary;
    end
end
p = 0.1;
n = 3;

a = zeros(1, 2^30);
b = zeros(1, 2^30);
c = zeros(1, 2^30);
d = zeros(1, 2^30);
% then OOM

sequences = generate_binary_sequences(n);

probs = zeros(1, 2^n);
for i = 1:2^n
    seq = sequences{i};
    ones_count = sum(seq == '1');
    probs(i) = p^ones_count * (1-p)^(n-ones_count);
end


% Calculate lengths
lengths = zeros(1, length(probs));
for i = 1:length(probs)
    lengths(i) = ceil(-log2(probs(i)));
end

lengths = sort(lengths);

max_length = max(lengths)
tree = create_tree(max_length + 1);
returns = 0;
% Store codes in a temporary cell array first
temp_codes = cell(1, length(lengths));
for i = 1:length(lengths)
    length_i = lengths(i);
    code = [];
    cur = 1; % MATLAB indices start at 1
    step = 0;
    while step < length_i
        fprintf("step %d at %d\n", step, cur);
        if tree(cur, 1) ~= -1
            if step == length_i - 1
                tree(cur, 1) = -1;
            end
            cur = tree(cur, 1);
            step = step + 1;
            code = [code, 0];
        elseif tree(cur, 2) ~= -1
            if step == length_i - 1
                tree(cur, 2) = -1;
                fprintf("deleting %d's connection from parent\n", cur)


                % go all the way to parent to waste
                % when we waste the second part we also waste the parent

                parent = ceil((cur-1)/2);
                if parent*2 == cur % left child (adjusting for 1-indexing)
                    tree(parent, 1) = -1;
                else % right child
                    tree(parent, 2) = -1;
                end

                delcur = cur;
                while delcur ~= 1
                    parent = ceil((delcur-1)/2);
                    if tree(parent) == -1 && tree(parent, 2) == -1
                        grandparent = ceil((parent-1)/2);
                        if grandparent*2 == parent % left child (adjusting for 1-indexing)
                            tree(grandparent, 1) = -1;
                        else % right child
                            tree(grandparent, 2) = -1;
                        end
                    else
                        break;
                    end
                    delcur = parent;
                end

            end
            cur = tree(cur, 2);
            step = step + 1;
            code = [code, 1];
        else
            % Find parent - adjusting for MATLAB's 1-indexed arrays
            parent = ceil((cur-1-1)/2) + 1;
            if parent*2 == cur % left child (adjusting for 1-indexing)
                tree(parent, 1) = -1;
            else % right child
                tree(parent, 2) = -1;
            end
            fprintf("returning\n");
            cur = 1;
            step = 0;
            code = [];
            returns = returns + 1;
        end
    end
    temp_codes{i} = code;
end

fprintf("had %d returns\n", returns);
% Convert to matrix format
code_lengths = cellfun(@length, temp_codes);
codes = zeros(length(temp_codes), max_length);
for i = 1:length(temp_codes)
    code = temp_codes{i};
    codes(i, 1:length(code)) = code;
end

length(codes)

[h, hl] = huffman_codebook(probs);

length(h)

for i = 1:size(codes, 1)
    code_str = mat2str(codes(i, 1:code_lengths(i)));
    fprintf('Code %d: %s (length: %d)\n', i, code_str, code_lengths(i));
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


function tree = create_tree(max_depth)
    % The "tree" is a matrix of nodes with columns representing left and right child indices
    % NULL is represented by -1
    
    % Calculate total number of nodes needed
    num_nodes = 2^max_depth - 1;
    tree = zeros(num_nodes, 2);
    
    % Create internal levels
    for i = 1:1:(2^(max_depth-1)-1)
        % Adjusting indices for MATLAB's 1-based indexing
        tree(i, :) = [i*2, i*2+1];
    end
    
    % Create bottom level with NULL children
    for i = 2^(max_depth-1):num_nodes
        tree(i, :) = [-1, -1];
    end
end


