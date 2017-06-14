%% Function to update the probability matrix

function P_mat = update_probability(probe_vec, Prob_matrix, result_vec)
Q = 1 - Prob_matrix;
tau = find(probe_vec); % Find the indices of the non-zero values in the probe vector
for i = 1: size(Prob_matrix, 1)
    w = prod(Q(i,tau)); % Calculate w_i values as in eq 12, section 2.5
% Update the probability values as in eq. 15, section 2.5
for j = tau
    if(result_vec(i) == 1)
        Prob_matrix(i,j) = Prob_matrix(i,j)/ (1 - w);
    else 
        Prob_matrix(i,j) = 0;
    end
end

end
P_mat = Prob_matrix;
end
