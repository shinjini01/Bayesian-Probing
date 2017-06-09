%% Function to find single best probes 
function tau = single_probe(prob_matrix, n)
Q = 1 - prob_matrix;
delta_U = max(sum(Q,1));
tau = find(sum(Q, 1) == delta_U, 1 ); %find minimum index with maximum value
a = 1:n; tau0 = tau;

   
while(length(a) ~= 0)
    a = find(1:n ~= tau0);
    for j = a
        b = (length(tau0)+1)*sum( prod(Q(:,tau0)).*Q(:,j))
        if(b >= delta_U)
            delta_U = b;
            tau = sort([tau, j]);
        end
    end
    if(tau == tau0)
        break;
    else
    tau0 = tau
    end
    
end

        
end


