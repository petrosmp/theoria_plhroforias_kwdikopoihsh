
function [codes_matrix, code_lengths] = shannon_fano_codebook_2(probs)
    lengths = zeros(1, length(probs));
    for i = 1:length(probs)
        lengths(i) = ceil(-log2(probs(i)));
    end
    
    lengths = sort(lengths);
    
    max_length = max(lengths);
    tree = create_tree(max_length + 1);
    
    temp_codes = cell(1, length(lengths));
    for i = 1:length(lengths)
        length_i = lengths(i);
        code = [];
        cur = 1;
        step = 0;
        while step < length_i
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
        
                    % when we waste the second part we also waste the parent
                    % go all the way to root to waste
                    parent = ceil((cur-1)/2);
                    if parent*2 == cur % left child
                        tree(parent, 1) = -1;
                    else % right child
                        tree(parent, 2) = -1;
                    end
    
                    delcur = cur;
                    while delcur ~= 1
                        parent = ceil((delcur-1)/2);
                        if tree(parent) == -1 && tree(parent, 2) == -1
                            grandparent = ceil((parent-1)/2);
                            if grandparent*2 == parent
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
                % this should never happen but we reset if it does
                cur = 1;
                step = 0;
                code = [];
                returns = returns + 1;
            end
        end
        temp_codes{i} = code;
    end
    
    code_lengths = cellfun(@length, temp_codes);
    codes_matrix = zeros(length(temp_codes), max_length);
    for i = 1:length(temp_codes)
        code = temp_codes{i};
        codes_matrix(i, 1:length(code)) = code;
    end
end

function tree = create_tree(max_depth)
    % The "tree" is a matrix of nodes with columns representing left and right child indices
    % NULL is represented by -1
    
    % Calculate total number of nodes needed
    num_nodes = 2^max_depth - 1;
    tree = zeros(num_nodes, 2);
    
    % Create internal levels
    for i = 1:1:(2^(max_depth)-1)
        % Adjusting indices for MATLAB's 1-based indexing
        tree(i, :) = [i*2, i*2+1];
    end
    
    % Create bottom level with NULL children
    for i = 2^(max_depth-1):num_nodes
        tree(i, :) = [-1, -1];
    end
end


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


