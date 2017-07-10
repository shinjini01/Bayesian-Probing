%% Function to find single best probe as in section 3.3
function tau = single_probe3(prob_matrix,m, n)
Q = 1 - prob_matrix;
delta_U = max(sum(Q,1)); % delta_U(P,tau) as defined in eq. 14, section 2.5
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index of columns with maximum value
indices = setdiff(1:n,tau); % Indices of the the columns not included in tau

while(~isempty(indices))
    indices = setdiff(1:n,tau);
    delta_U_i = zeros(1,length(indices));
    indices_U_i = cell(1, length(indices));
    
    
    for j = 1:length(indices)
        columns = nchoosek(indices, j); %matrix of column index combinations
        delta_U_k = 0;indices_i0 = 0;
        for k = 1:size(columns, 1)
            indices_k = columns(k,:);
            
            delta_U_i_k = 0;tau_i = 0;
            % recalculating delta_U for additional columns
            for i = 1:m
                
                indices_i = indices_k(Q(i, indices_k) <1); %column indices such that p_ij >0 in the ith row
                tau_i = tau(Q(i, tau) <1);
                
                
                delta_U_i_k = delta_U_i_k + (length(tau_i)+length(indices_i))*sum( prod(Q(i,tau_i),2)*prod(Q(i,indices_i),2));
            end
                if(delta_U_i_k > delta_U_k)
                    delta_U_k = delta_U_i_k;
                    indices_i0 = unique(sort(horzcat(tau_i, columns(k,:))));
                end
            end
            delta_U_i(j) = delta_U_k;
            indices_U_i{j} = indices_i0;
        end
        
    % extending/modifying tau, i.e., support of the probe vector on the
    % basis of changes in delta_U
    if(max(delta_U_i) > delta_U)
        delta_U = max(delta_U_i);
        max_index = find(delta_U_i == max(delta_U_i),1);
        tau = unique(sort(horzcat(tau, indices_U_i{max_index})));
    else
        break
    end
    
end


end
