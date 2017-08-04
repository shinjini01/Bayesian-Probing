
function n_itr2 = find_probe(P,S)
[m_s, n_s] = size(S);
n_itr2 = 0;
P_prime = zeros(m_s, n_s); S_prime =S;
tic
while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S, S<1))))
    [m n ] = size(P); indices = NaN(1,n);
    probe_mat = zeros(n, 1);r_mat = zeros(m, 1);
    
    tau = single_probe4(P,n);
    probe_mat(tau(1:length(tau)),:) = 1
    
    
    % The corresponding result matrix
    r_mat(find(S_prime* probe_mat)) = 1;
    
     P = update_probability(probe_mat, P, r_mat)
%     for i = 1:n
%         for j = 1:n_s
%             if(sum(abs(P(:,i)-S(:,j))< 0.0000001) == m)
%                 indices(i) = i;
%                 P_prime(:,j) = P(:,i);
%             end
%         end
%     end
%     
%     %round(P_prime(80:100, 80:100))
%     % Remove the fully identified columns from the probability matrix
%     % before the next iteration
%     col_indices = (1:n).*(isnan(indices));
%     P = P(:,col_indices(col_indices>0));
%     % Remove the corresponding columns from the sparsity matrix
%     S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr2 = n_itr2+1
    
end
toc
end
