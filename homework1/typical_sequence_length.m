EPSILON = 0.05;
p = 0.45;

binary_source_entropy = @(p) -p*log2(p) - (1-p)*log2(1-p);

H = binary_source_entropy(p);

A_ns = [];
upper_bounds = [];
lower_bounds = [];
n_values = 10:10:1000;

for n = n_values
    A_n = nchoosek(n, round(n*p));
    
    upper_bound = 2^(n*(H + EPSILON));
    lower_bound = (1 - EPSILON) * 2^(n*(H - EPSILON));
    
    A_ns = [A_ns, A_n];
    upper_bounds = [upper_bounds, upper_bound];
    lower_bounds = [lower_bounds, lower_bound];
    
    disp(n)
end


figure;
hold on;
ax=gca;
set(ax, 'YScale', 'log');

semilogy(n_values, A_ns, 'b', DisplayName="A_n (p=0.45)");
semilogy(n_values, upper_bounds, 'r', DisplayName="Upper Bound (p=0.45)");
semilogy(n_values, lower_bounds, 'g', DisplayName="Lower Bound (p=0.45)");
hold off;

xlabel('n');
ylabel('A_n, Upper Bound, Lower Bound');
title('A_n, Upper Bound, and Lower Bound');
grid on;




p = 0.1;
H = binary_source_entropy(p);

A_ns = [];
upper_bounds = [];
lower_bounds = [];
n_values = 10:10:1000;


for n = n_values
    A_n = nchoosek(n, round(n*p));
    
    upper_bound = 2^(n*(H + EPSILON));
    lower_bound = (1 - EPSILON) * 2^(n*(H - EPSILON));
    
    A_ns = [A_ns, A_n];
    upper_bounds = [upper_bounds, upper_bound];
    lower_bounds = [lower_bounds, lower_bound];
    
    disp(n)
end


hold on;
ax=gca;
set(ax, 'YScale', 'log');

semilogy(n_values, A_ns, 'black', DisplayName="A_n (p=0.1)");
semilogy(n_values, upper_bounds, 'm', DisplayName="Upper Bound (p=0.1)");
semilogy(n_values, lower_bounds, 'c', DisplayName="Lower Bound (p=0.1)");
hold off;

xlabel('n');
ylabel('A_n, Upper Bound, Lower Bound');
legend();
title('A_n, Upper Bound, and Lower Bound');
grid on;
