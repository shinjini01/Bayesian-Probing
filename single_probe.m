%% Function to find single best probes 
function tau = single_probe(prob_matrix, n)
Q = 1 - prob_matrix;
delta_U = max(sum(Q,1));
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index with maximum value
indices = find(1:n ~= tau);

   
while(~isempty(indices))
    w_i = [];
    for j = 1:length(indices)
        w_i = [w_i,(length(tau)+j)*sum( prod(Q(:,tau),2).*Q(:,indices[1:j]))];
        end;
        
        if(max(w_i) >= delta_U)
            delta_U = w_i;
            tau = sort([tau, 1:find(w_i==max(w_i))]);
        end
    else
    break;
    end
    
end

        
end


