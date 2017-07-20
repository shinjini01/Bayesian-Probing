%% Function to update the probability matrix
%  probe_mat is a matrix where each column comprises a probe. Order: n x K
%  result_mat is a matrix of the r values where each row corresponds to a
%  row of S multiplied with each of the probes in the probe_mat. 
%  Order: m x K
%  Prob_matrix is the matrix of probabilities

function P_mat = update_probability_bundle(probe_mat, P, r_mat)
Q = 1 - P; tau = cell(1, size(probe_mat,2));
for i = 1:size(probe_mat,2)
    tau{i} = find(probe_mat(:,i)); % Find the indices of the non-zero values in the probe vectors
end


for i = 1: (size(P, 1)) %Taking one row of probability matrix
    result_vec = r_mat(i,:);
    for j = 1:size(probe_mat,2)
        w(j) = prod(Q(i,tau{j})); % Calculate w_i values as in eq 12, section 2.5
        
        
        %     result_0 = find(result_vec == 0);
        %     result_1 = find(result_vec == 1);
        
        % Update the probability values as in eq. 15, section 2.5
        for k = tau{j}
            if(result_vec(j) == 1)
               if(P(i,k)>1)
                   P(i,k) = P(i,k)/P(i,k);
                  
                else 
                   P(i,k) = P(i,k)/ (1 - w(j));
               end
            else
                P(i,k) = 0;
            end
            
        end
    end
end

P_mat = P;
end