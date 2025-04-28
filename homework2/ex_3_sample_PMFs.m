n=10;
p = 0.1;

fileID = fopen("original_data_5.bin", 'rb');
X = fread(fileID, 'ubit1');
fclose(fileID);

fileID = fopen("H_compressed_data_5_10.bin", 'rb');
Y_h = fread(fileID, 'ubit1');
fclose(fileID);

fileID = fopen("SF_compressed_data_5_10.bin", 'rb');
Y_sf = fread(fileID, 'ubit1');
fclose(fileID);

X_chunks = floor(length(X) / n);
Y_h_chunks = floor(length(Y_h) / n);
Y_sf_chunks = floor(length(Y_sf) / n);

X_reshaped = reshape(X(1:X_chunks * n), n, []).';
Y_h_reshaped = reshape(Y_h(1:Y_h_chunks * n), n, []).';
Y_sf_reshaped = reshape(Y_sf(1:Y_sf_chunks * n), n, []).';

X_sequences_decimal = bi2de(X_reshaped, 'left-msb');
Y_h_sequences_decimal = bi2de(Y_h_reshaped, 'left-msb');
Y_sf_sequences_decimal = bi2de(Y_sf_reshaped, 'left-msb');

% Count the occurrences of each possible sequence (0 to 1023)
X_pmf_counts = histcounts(X_sequences_decimal, 0:1024);
X_pmf = pmf_counts / sum(X_pmf_counts); % Normalize to get PMF

Y_h_pmf_counts = histcounts(Y_h_sequences_decimal, 0:1024);
Y_h_pmf = Y_h_pmf_counts / sum(Y_h_pmf_counts); % Normalize to get PMF

Y_sf_pmf_counts = histcounts(Y_sf_sequences_decimal, 0:1024);
Y_sf_pmf = Y_sf_pmf_counts / sum(Y_sf_pmf_counts); % Normalize to get PMF

function res = H(p)
    res = -p * log2(p) - (1-p) * log2(1-p);
end

figure;
semilogy(0:1023, sort(X_pmf), DisplayName="Original", LineWidth=1.5);
hold on;
semilogy(0:1023, sort(Y_h_pmf), DisplayName="Huffman", LineWidth=1.5);
semilogy(0:1023, sort(Y_sf_pmf), DisplayName="Shannon-Fano", LineWidth=1.5);
legend();
xlabel('10-bit sequence (decimal representation)');
ylabel('Probability');
title('PMF of 10-bit Sequences');
grid on;
xlim([0 1023]);

% EPSILON = 0.51;
% typical_sequence_lower = 2^(-n*(H(p)+EPSILON));
% typical_sequence_upper = 2^(-n*(H(p)-EPSILON));
% yline(typical_sequence_lower,  HandleVisibility="off");
% yline(typical_sequence_upper, HandleVisibility="off");
