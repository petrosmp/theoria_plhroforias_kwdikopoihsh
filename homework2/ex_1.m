%{
Source with alphabet X={1,2,3,4} with probabilities 1/2, 1/4, 1/8, 1/8.

Huffman code (derived on paper):
    1 -> 0
    2 -> 10
    3 -> 110
    4 -> 111

The script runs a number of experiments in each of which a sequence of
length n is produced from the source, coded with the Huffman code above,
decoded and compared with the initial. Statistics are also gathered
regarding average code word length and code word distributions.
%}

clear all
close all

X = [1, 2, 3, 4];
source_probs = [1/2, 1/4, 1/8, 1/8];
entropy = 1.75;
n = 100;
codebook = {'0', '10', '110', '111'};
num_of_experiments = 1000;


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

errors = 0;
avg_code_word_lengths = zeros(1, num_of_experiments);
d_5 = cell(1, num_of_experiments);
d_6 = cell(1, num_of_experiments);
joint = cell(1, num_of_experiments);

for exp = 1:num_of_experiments
    raw = randsample(X, n, true, source_probs); % needs statistics toolbox
    
    encoded = huffman(raw, codebook);
    decoded = dehuffman(encoded, codebook);
    
    avg_code_word_lengths(exp) = length(encoded) / n;
    
    d_5{exp} = encoded(5);
    d_6{exp} = encoded(6);
    joint{exp} = [encoded(5), encoded(6)];

    if length(raw) ~= length(decoded)
        fprintf('Error: Length mismatch in experiment %d\n', exp);
        errors = errors + 1;
    else
        for i = 1:length(raw)
            if raw(i) ~= decoded(i)
                fprintf('Error at position %d, raw [%d] != decoded [%d]\n', i, raw(i), decoded(i));
                errors = errors + 1;
            end
        end
    end
end
fprintf('Total errors: %d\n', errors);


% aggregate results for plotting
unique_lengths = unique(avg_code_word_lengths);
length_counts = histcounts(avg_code_word_lengths, [unique_lengths, max(unique_lengths)+0.0001]);
length_counts = length_counts / num_of_experiments;

d_5_counts = containers.Map({'0', '1'}, [0, 0]);
d_6_counts = containers.Map({'0', '1'}, [0, 0]);
joint_counts = containers.Map({'00', '01', '10', '11'}, [0, 0, 0, 0]);

for i = 1:num_of_experiments
    if ~isempty(d_5{i})
        d_5_counts(d_5{i}) = d_5_counts(d_5{i}) + 1;
    end
    
    if ~isempty(d_6{i})
        d_6_counts(d_6{i}) = d_6_counts(d_6{i}) + 1;
    end
    
    if ~isempty(joint{i})
        if isKey(joint_counts, joint{i})
            joint_counts(joint{i}) = joint_counts(joint{i}) + 1;
        end
    end
end

keys_d5 = keys(d_5_counts);
keys_d6 = keys(d_6_counts);
keys_joint = keys(joint_counts);


% plot everything

figure;
plot(unique_lengths, length_counts);
title('Distribution of average code word length');
hold on;
xline(entropy, '--k', 'entropy', LabelOrientation='horizontal');

figure;
subplot(2, 2, 1);
bar(1:2, [d_5_counts('0'), d_5_counts('1')]);
title('Distribution of d_5');
xticks(1:2);
xticklabels({'0', '1'});

subplot(2, 2, 2);
bar(1:2, [d_6_counts('0'), d_6_counts('1')], FaceColor=[1, 0.5, 0]);
title('Distribution of d_6');
xticks(1:2);
xticklabels({'0', '1'});

subplot(2, 1, 2);
joint_values = zeros(1, 4);
for i = 1:length(keys_joint)
    if strcmp(keys_joint{i}, '00')
        joint_values(1) = joint_counts(keys_joint{i});
    elseif strcmp(keys_joint{i}, '01')
        joint_values(2) = joint_counts(keys_joint{i});
    elseif strcmp(keys_joint{i}, '10')
        joint_values(3) = joint_counts(keys_joint{i});
    elseif strcmp(keys_joint{i}, '11')
        joint_values(4) = joint_counts(keys_joint{i});
    end
end
bar(1:4, joint_values, 'FaceColor', 'g');
title('Joint distribution of d_5 and d_6');
xticks(1:4);
xticklabels({'00', '01', '10', '11'});
