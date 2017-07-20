function probe_matrix = cluster_probe(P,m, n, K)
Q = 1-P; probe_matrix = zeros(n, K); indices = 1:n;
col_sums = zeros(1,n);
for j = 1:n
    unit_indices = find(Q(:,j) < 1);
    col_sums(j) = sum(Q(unit_indices, j));
end

clusters = kmeans(transpose(col_sums), K, 'MaxIter', 100000);
for i = 1:K
    tau = clusters == i;
    probe_matrix(tau,i ) = 1;
end


tau = zeros(n/K, K);
[B, I] =  sort(col_sums, 'descend');
tau(1,:) = I(1:K); % First element in each group selected

indices(tau(1,:)) = [];


for i = 1:(n/K - 1)
    for j = 1:K
        delta_U_j = 0;
        for k = 1:length(indices)
            mat = Q(:,[tau((1:i),j) indices(k)]);
            delta_U = sum( (sum((mat(:,tau((1:i),j)))<1, 2) + (mat(:,indices(k)) <1)) .* ( prod(mat(:,tau((1:i),j)),2) .* mat(:,indices(k)) )) ;
            if( delta_U > delta_U_j )
                tau((i+1),j) = k;
                delta_U_j = delta_U;
            end
        end
    end
end

probe_matrix        



end



s(:, [t(1,(1:3)) 5])
(sum(s(:,t(1,(1:3)))<1,2) + (s(:,5) <1)) .* ( prod(s(:,t(1,(1:3))),2) .* s(:,5))