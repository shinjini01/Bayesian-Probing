%% Function to find single best probe as in section 3.3
function tau = single_probe4(P)
Q = 1 - P ;
sum_indices = sum( (Q<1) .* Q, 1);
delta_U = max(sum_indices);  %delta_U0 = delta_U;% delta_U(P,tau) as defined in eq. 14, section 2.5
tau = find(sum_indices == delta_U,1);
[~, I] =  sort(sum_indices, 'descend');
indices = setdiff(I,tau); % Indices of the the columns not included in tau

while(~isempty(indices))
    
    indices = setdiff(I,tau);
    %matrix of column index combinations
    delta_U_k = 0; indices_i0 = 0;
    l_tau = sum(Q(:,tau)<1, 2); w_tau = prod(Q(:,tau),2);
    
    for k = 1:length(indices)
        
        % recalculating delta_U for additional columns
        delta_U_i_k = sum( (l_tau + sum(Q(:,indices(k)) <1,2)) .* ( w_tau .* prod(Q(:,indices(k)),2) ));
        
        if(delta_U_i_k > delta_U_k)
            delta_U_k = delta_U_i_k ; %Maximum delta value from adding a column to a probe
            indices_i0 = unique(sort(horzcat(tau, indices(k))));
        end
    end
    
    if(delta_U_k > delta_U)
        delta_U = delta_U_k;
        tau = indices_i0; % determining the support of the probe vector
    else
        break
    end    
end

end
