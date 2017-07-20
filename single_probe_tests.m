%% Bundle Probing %%
load can_96.mat
sparsemat = full(Problem.A);
S = abs(sign(sparsemat));
S = diag(ones(10,1));    % Diagonal structure
S(:,100) = 1; S(100,:) = 1;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 1; batch_size = round(n_s/K);
%P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
P = rand(m_s, n_s);
mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
P = min(mnzi,nnzj)/n_s;

P = random('beta', 0.5, (1 - (1/(2*log(n_s))))*0.5, n_s);

P_prime2 = zeros(m_s, n_s);
P = ones(m_s, n_s)/n_s;
n_itr2 = zeros(1,5); %probe_matrices = cell(5,1);
tic

for l = 1:5
   S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 1; batch_size = round(n_s/K); 
if(l ==1)
    P = ones(m_s, n_s)/n_s;
elseif (l==2)
    P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
elseif (l==3)
    P = rand(1/n_s, (1-1/n_s), n_s);
elseif (l==4)
    mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
    P = min(mnzi,nnzj)/n_s;
elseif (l==5)
    a = random('beta', (1 + 1/n_s), (1-1/n_s), m_s, n_s, 500);
    P = mean(a, 3);
end
P_prime2 = zeros(m_s, n_s);  n_itr2 = zeros(1,5);  
while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
    
    P_prime = P;
    [m_p, n_p] = size(P);
    probe_mat = zeros(n_p, K);r_mat = zeros(m_p, K);
    indices = NaN(1,n_p);ind = 1:n_p;
    
    
    % Finding the bundle of probes
    
    %      [m,n] = size(P_prime);
    %         tau = single_probe3(P_prime,m,n);
    %         probe_mat(tau) = 1;
    %
    %         P_prime(:,tau(1:size(tau,2)))= ones(m,size(tau,2));
    
    for i = 1:K
        [m,n] = size(P_prime);
        tau = single_probe4(P_prime,m,n);
        probe_mat(tau(1:min(size(tau,2),batch_size)),i) = 1;
        
         P_prime(:,tau(1:min(size(tau,2),batch_size)))= ones(m,min(size(tau,2),batch_size));
        
        
    end
    % The corresponding result matrix
    r_mat(find(S_prime * probe_mat)) = 1;
    
    
    if(K>1)
        P = update_probability_bundle(probe_mat, P, r_mat);
    else
        P = update_probability(probe_mat, P, r_mat);
    end
    
    
    for i = 1:n_p
        for j = 1:n_s
            if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == m_p)
                indices(i) = i;
                P_prime2(:,j) = P(:,i);
            end
        end
    end
    % round(P_prime2(480:500,480:500))
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr2= n_itr2+1;
    probe_matrices2{n_itr2}= probe_mat;
end
toc
cellfun(@sum, probe_matrices)


%% Second Bundle Probing method implemented
S = diag(ones(100,1));    % Diagonal structure
S(:,100) = 1; S(100,:) = 1;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 10; batch_size = round(n_s/K);
%P = diag(sum(S,2)/n_s)*ones(m_s, n_s);

mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
P = min(mnzi,nnzj)/n_s;

P = random('beta', 0.5, (1 - (1/(2*log(n_s))))*0.5, n_s);
P = ones(m_s, n_s)/n_s;
P_prime2 = zeros(m_s, n_s);

n_itr2 = 0; probe_matrices2 = cell(2,1);

tic
while( sum(sum(and(0 < P, P <1))) ~= sum(sum(and(0<S_prime, S_prime<1))))
    
    P_prime = P;
    [m_p, n_p] = size(P);
    r_mat = zeros(m_p, K);
    indices = NaN(1,n_p);ind = 1:n_p;
    probe_mat = bundle_probe2(P, n_p, K);
    
   
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
    % round(P_prime2(80:96,80:96))
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr2= n_itr2+1;
    probe_matrices{n_itr2}= probe_mat;
end
toc
cellfun(@sum, probe_matrices)