
function n_itr2 = find_probe(P,S)
[m_s, n_s] = size(S);
P_prime2 = zeros(m_s, n_s);
n_itr2 = 0;
S_prime = S;

while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
    
    P_prime = P;
    [m_p, n_p] = size(P);
    probe_mat = zeros(n_p, 1);r_mat = zeros(m_p, 1);
    indices = NaN(1,n_p);ind = 1:n_p;
    
    
    
    tau = single_probe4(P_prime,m_p,n_p);
    probe_mat(tau(1:size(tau,2)),:) = 1;
    
    
    % The corresponding result matrix
    r_mat(find(S_prime * probe_mat)) = 1;
    
    P = update_probability(probe_mat, P, r_mat);
    
    
    for i = 1:n_p
        for j = 1:n_s
            if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == m_p)
                indices(i) = i;
                P_prime2(:,j) = P(:,i);
            end
        end
    end
    % round(P_prime2(80:96,80:96))
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    %n_itr2(l)= n_itr2(l)+1;
    n_itr2 = n_itr2+1;
    
end
end
