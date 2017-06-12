%% Function to find single best probes 
function tau = single_probe2(prob_matrix, n)
Q = 1 - prob_matrix; 
delta_U = max(sum(Q,1));
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index with maximum value
indices = setdiff(1:n,tau);
   
while(~isempty(indices))
    indices = setdiff(1:n,tau);
    delta_U_i = NaN(1,length(indices));
    for j = 1:length(indices)
        delta_U_i(j) = (length(tau)+j)*sum( prod(Q(:,tau),2).*prod(Q(:,indices(1:j)),2));
    end
    delta_U_i = delta_U_i(~isnan(delta_U_i));
        if(max(delta_U_i) > delta_U)
            delta_U = max(delta_U_i); max_index = find(delta_U_i == max(delta_U_i),1);
            tau = unique(sort(horzcat(tau, indices(1:max_index))));
        else
            break
        end
      
end

        
end

