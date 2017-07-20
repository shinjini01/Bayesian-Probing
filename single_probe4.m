%% Function to find single best probe as in section 3.3
function tau = single_probe4(P,m, n)
Q = 1 - P; 
sum_indices = sum( (Q<1) .* Q, 1);
delta_U = max(sum_indices); % delta_U(P,tau) as defined in eq. 14, section 2.5
tau = find(sum_indices == delta_U,1);
indices = setdiff(1:n,tau); % Indices of the the columns not included in tau

while(~isempty(indices))
    indices = setdiff(1:n,tau);
     
    
    columns = nchoosek(indices, 1); %matrix of column index combinations
    delta_U_k = 0; indices_i0 = 0;
    for k = 1:size(columns, 1)
             
        % recalculating delta_U for additional columns
        delta_U_i_k = sum( (sum(Q(:,tau)<1, 2) + sum(Q(:,columns(k,:)) <1,2)) .* ( prod(Q(:,tau),2) .* prod(Q(:,columns(k,:)),2) ));
        
        if(delta_U_i_k > delta_U_k)
            delta_U_k = delta_U_i_k;  %Maximum delta value from a
            indices_i0 = unique(sort(horzcat(tau, columns(k,:))));
            
        end
    end
   
    if(delta_U_k > delta_U)
        delta_U = delta_U_k;
      
        tau = unique(sort(horzcat(tau, indices_i0)));
    else
        break
    end
    
end


end

