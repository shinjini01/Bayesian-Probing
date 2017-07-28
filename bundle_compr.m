%% Bundle Probing %%
% load can_96.mat
% sparsemat = full(Problem.A);
% S = abs(sign(sparsemat));
 %S = diag(ones(100,1));    % Diagonal structure
 S = speye(1000);
% S(:,100) = 1; S(100,:) = 1;
K = [2 4 16 ];  n_itr2 = zeros(1,length(K));
K = 32; n_itr2 = 0;


    
    probe_matrices2 = cell(2,1)
    
    for l = 1:length(K)
    S_prime = S;
    [m_s, n_s] = size(S);     % Size of the sparse matrix
    batch_size = round(n_s/K(l));
    P_prime2 = zeros(m_s, n_s);
    P = ones(m_s, n_s)/n_s;
    tic
    while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
        
        P_prime = P;
        [m_p, n_p] = size(P);
        probe_mat = zeros(n_p, K(l));r_mat = zeros(m_p, K(l));
        indices = NaN(1,n_p);ind = 1:n_p;
        
        
        for i = 1:K(l)
            [m,n] = size(P_prime);
            tau = single_probe4(P_prime,n);
            probe_mat(tau(1:min(size(tau,2),batch_size)),i) = 1;
            
            P_prime(:,tau(1:min(size(tau,2),batch_size)))= ones(m,min(size(tau,2),batch_size));
            
            
        end
        % The corresponding result matrix
        r_mat(find(S_prime * probe_mat)) = 1;
        
        
        
        P = update_probability_bundle(probe_mat, P, r_mat);
        
        
        for i = 1:n_p
            for j = 1:n_s
                if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == m_p)
                    indices(i) = i;
                    P_prime2(:,j) = P(:,i);
                end
            end
        end
        % round(P_prime2(980:1000,980:1000))
        
        % Remove the fully identified columns from the probability matrix
        % before the next iteration
        col_indices = (1:n_p).*(isnan(indices));
        P = P(:,col_indices(col_indices>0));
        % Remove the corresponding columns from the sparsity matrix
        S_prime = S_prime(:,col_indices(col_indices>0));
        
        n_itr2(l)= n_itr2(l)+1;
        %probe_matrices2{n_itr2}= probe_mat;
    end
    toc
    n_itr2(l)
    end
    for h = 1:n_itr2
        a = sum(probe_matrices2{h},1); sum(a~=0)
    end





%% Second Bundle Probing method implemented
% S = diag(ones(100,1));    % Diagonal structure
% S(:,100) = 1; S(100,:) = 1;K  = 8
 n_itr3 = zeros(1, length(K));

    for l = 1:length(K)
    S_prime = S;
    [m_s, n_s] = size(S);     % Size of the sparse matrix
    %K = 2; batch_size = round(n_s/K);
    %P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
    
    %mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
    %P = min(mnzi,nnzj)/n_s;
    
    %P = random('beta', 0.5, (1 - (1/(2*log(n_s))))*0.5, n_s);
    P = ones(m_s, n_s)/n_s;
    P_prime2 = zeros(m_s, n_s);
    
    %probe_matrices3 = cell(2,1);
    
    tic
    while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
             
         r_mat = zeros(size(P,1), K(l));
        indices = NaN(1,size(P,2));
        probe_mat = bundle_probe2(P, size(P,2), K(l));
        
        
        % The corresponding result matrix
        r_mat(find(S_prime * probe_mat)) = 1;
        
        
        P = update_probability_bundle2(probe_mat, P, r_mat);
        
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
        n_itr3(l)= n_itr3(l)+1;
       
    end
    toc
    n_itr3(l)
    end
    x = ['n_itr3 = ', num2str(n_itr3(p))]; disp(x)
    for h = 1:n_itr3(p)
        a = sum(probe_matrices3{h},1); sum(a~=0)
    end

%% Third Bundle probing method

P1 = ones(m_s, n_s)/min(m_s, n_s); P2 = ones(m_s, n_s)/2; P3 = 1-P1;
S_prime = S;
[m_s, n_s] = size(S); K = 10; batch_size = round(n_s/K);