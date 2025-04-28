SF_sizes = zeros(10, 1);
H_sizes = zeros(10, 1);
upper_bound = zeros(10, 1);
p = 0.1;
fileID = fopen("original_data_5.bin", 'rb');
A = fread(fileID, 'ubit1');
M = length(A);
fclose(fileID);

function res = H(p)
    res = -p * log2(p) - (1-p) * log2(1-p);
end


for n=1:10

    SF_filename = sprintf('SF_compressed_data_5_%dBINEXP.bin', n);
    H_filename = sprintf('H_compressed_data_5_%d.bin', n);

    upper_bound(n) = (H(p) + 1/n) * M;

    if exist(SF_filename, 'file')
        fileID = fopen(SF_filename, 'rb');
        A = fread(fileID, 'ubit1');
        fclose(fileID);
        SF_sizes(n) = length(A);
    else
        error("Error! File '%s' does not exist\n", SF_filename);
    end
    
    if exist(H_filename, 'file')
        fileID = fopen(H_filename, 'rb');
        A = fread(fileID, 'ubit1');
        fclose(fileID);
        H_sizes(n) = length(A);
    else
        error("Error! File '%s' does not exist\n", H_filename);
    end
    
end

x = 1:10;

figure;
plot(x, SF_sizes, DisplayName='Shannon-Fano', LineWidth=1.5);
hold on;
plot(x, H_sizes, DisplayName='Huffman', LineWidth=1.5);
yline(H(p)*M, DisplayName='Lower Bound', LineWidth=1.5);
plot(x, upper_bound, '--', DisplayName='Upper Bound', Color="black", LineWidth=1);
legend();

xlabel('n');
ylabel('compressed file size (in bits)');
title('Comparison of SF and H Compressed File Sizes');
