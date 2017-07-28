% Function to bundle probing in a new way. The largest column sums are
% allotted to different K batches. The remaining columns are searched to
% find the index that increases delta_U for the first batch. The remaining
% columns are again searched through to find the best second candidate for
% the second batch and so on, till all column  indices are exhausted and
% the column indices are sorted into a batch_size x K matrix.

function probe_matrix = bundle_probe2(P,n, K)
Q = 1-P; probe_matrix = zeros(n, K); indices = 1:n;
col_sums = sum( (Q<1) .* Q, 1); batch_size = round(n/K); 


tau = zeros(batch_size, K);
[~, I] =  sort(col_sums, 'descend');
tau(1,(1:min(n,K))) = I(1:min(n,K)); % First element in each group selected

indices(nonzeros(tau(1,:))) = [];    % Column indices not included in the first set 
delta_U = col_sums(nonzeros(tau(1,:)));

for i = 1:(batch_size - 1)
    for j = 1:K
        k_star = 0;
        for k = 1:length(indices)
            
            delta_U_j = sum( (sum((Q(:,nonzeros(tau((1:i),j))))<1, 2) + (Q(:,indices(k)) <1)) .* ( prod(Q(:,nonzeros(tau((1:i),j))),2) .* Q(:,indices(k)) ));
            if( delta_U_j > delta_U(j) )
                tau((i+1),j) = indices(k);
                k_star = k;
                delta_U(j) = delta_U_j;
            end
        end
        indices(nonzeros(k_star)) = [];
    end
        
    
end

for j = 1:K
    probe_matrix(nonzeros(tau(:,j)),j) = 1;
end



end
