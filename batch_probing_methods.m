
% The two bundle probing methods, i.e., Griewank-Mitev's method and 
% the proposed method can be executed using the following code
S = diag(ones(8,1));    % Diagonal structure
%S = speye(1000);

K = 2; n_itr = 0;      % K: # of probes in a batch; n_itr: # of batches

S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
probe_size = round(n_s/K);
P_prime2 = zeros(m_s, n_s);
P = ones(m_s, n_s)/n_s;   % Initializing the probability matrix

while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
    
    P_prime = P;
    [m_p, n_p] = size(P);
    probe_mat = zeros(n_p, K);r_mat = zeros(m_p, K);
    indices = NaN(1,n_p);    
    
    for i = 1:K
        [m,n] = size(P_prime);
        tau = single_probe4(P_prime,n);
        probe_mat(tau(1:min(size(tau,2),probe_size)),i) = 1;        
        P_prime(:,tau(1:min(size(tau,2),probe_size)))= ones(m,min(size(tau,2),probe_size));        
        
    end
    % The corresponding result matrix
    r_mat(find(S_prime * probe_mat)) = 1;
    
    % Updating the probability matrix using r_mat and probe_mat    
    P = update_probability_bundle(probe_mat, P, r_mat);
    
    
    % Collecting indexes of the already discovered columns in order to 
    % remvove them from further consideration and quicken the process
    for i = 1:n_p
        for j = 1:n_s
            if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == m_p)
                indices(i) = i;
                P_prime2(:,j) = P(:,i);
            end
        end
    end  
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr= n_itr+1;
    
end

n_itr;
P_prime2;





%% Second Bundle Probing method implemented
n_itr2 = 0;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix

probe_size = round(n_s/K);

P = ones(m_s, n_s)/n_s;
P_prime2 = zeros(m_s, n_s);

while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
    
    r_mat = zeros(size(P,1), K);
    indices = NaN(1,size(P,2));
    probe_mat = bundle_probe2(P, size(P,2), K);    
    
    % The corresponding result matrix
    r_mat(find(S_prime * probe_mat)) = 1;    
    
    % Updating the probability matrix using r_mat and probe_mat    
    P = update_probability_bundle(probe_mat, P, r_mat);
    
    % Collecting indexes of the already discovered columns in order to 
    % remvove them from further consideration and quicken the process
    for i = 1:size(P,2)
        for j = 1:n_s
            if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == size(P,1))
                indices(i) = i;
                P_prime2(:,j) = P(:,i);
            end
        end
    end
    
    col_indices = (1:size(P,2)).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    
    S_prime = S_prime(:,col_indices(col_indices>0));
    n_itr2= n_itr2+1;
    
    
end

n_itr2;
P_prime2;   
