function tau = single_probe(prob_matrix, n)
Q = 1 - prob_matrix; 
delta_U = max(sum(Q,1));
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index with maximum value
indices = setdiff(1:n,tau);
   
while(~isempty(indices))
    indices = setdiff(1:n,tau);
    w_i = NaN(1,length(indices));
    for j = 1:length(indices)
        w_i(j) = (length(tau)+j)*sum( prod(Q(:,tau),2).*prod(Q(:,indices(1:j)),2));
    end
    w_i = w_i(~isnan(w_i));
        if(max(w_i) > delta_U)
            delta_U = max(w_i); max_index = find(w_i == max(w_i),1);
            tau = unique(sort(horzcat(tau, indices(1:max_index))));
        else
            break
        end
      
end

end
