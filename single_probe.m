function tau = single_probe(prob_matrix,m, n)
Q = 1 - prob_matrix; 
delta_U = max(sum(Q,1));
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index with maximum value
indices = setdiff(1:n,tau);
   
while(~isempty(indices))
    indices = setdiff(1:n,tau);
    delta_U_i = zeros(1,length(indices));
   
    for j = 1:length(indices)
        
        for i = 1:m
            indices_i = indices(Q(i,indices(1:j))<1);
            tau_i = tau(Q(i, tau) <1);
            delta_U_i(j) = delta_U_i(j) + (length(tau_i)+length(indices_i))*sum( prod(Q(i,tau_i),2)*prod(Q(i,indices_i),2));
        end
    end
    %delta_U_i = delta_U_i(~isnan(delta_U_i));
        if(max(delta_U_i) > delta_U)
            delta_U = max(delta_U_i); max_index = find(delta_U_i == max(delta_U_i),1);
            tau = unique(sort(horzcat(tau, indices(1:max_index))));
        else
            break
        end
      
end

        
end

