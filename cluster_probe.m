function probe_matrix = cluster_probe(P, n, K)
Q = 1-P; probe_matrix = zeros(n, K);
col_sums = sum(Q, 1);
clusters = kmeans(transpose(col_sums), K, 'MaxIter', 100000, 'Display', 'iter');
for i = 1:K
    tau = find(clusters == i);
    probe_matrix(tau,i ) = 1;
end
end



