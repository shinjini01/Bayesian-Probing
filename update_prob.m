% P = ones(8,8)/8 ;    % Initialize the probability matrix
% Q = 1 - P ;          % 1 - Probability matrix, required to calculate w values
% % To initialize the matrix with each column as a probe vector
% a1 = [ones(4,1); zeros(4,1)];
% a2 = repmat([1;1;0;0],2,1);
% a3 = repmat([1;0],4,1);
% t = [a1, 1-a1, a2, 1 - a2, a3, 1 - a3]; % t is the required matrix of probes
% a0 = repmat([1;0],3,1);
% a4 = [a0; 1;1];
% a5 = [1-a0; 1;1];
% r.mat = [a1, 1 - a1, a2, 1 - a2, a4, a5]; % r.mat is the result matrix
%
%
%
% for k = 1:size(t,2)
% P = update_probability(t(:,k), P, r.mat(:,k))
%
% end
%
%%%% Generalize the Algorithm %%%%

%    Define the sparsity matrix

S = diag(ones(64,1));    % Diagonal structure
S(:,64) = 1; S(64,:) = 1; % Arrow type structure

%%% Sparsity matrix from the suitesparse.com website
%   Weblink: http://yifanhu.net/GALLERY/GRAPHS/search.html

%load cage4.mat
load can_96.mat
sparsemat = full(Problem.A);
S = abs(sign(sparsemat)); % The sparse matrix

[m_s, n_s] = size(S);     % Size of the sparse matrix
S_prime = S;
%P = ones(m_s,n_s)/min(m_s,n_s); % We begin with P matrix with identical elements 1/min(m, n)
%P  = ones(m_s, n_s)/2
P = diag(sum(S,2)/n_s)*ones(m_s, n_s)
P_prime = zeros(m_s,n_s);
n_itr1 =0;                 % Counter for the number of interations

tic
while( sum(and(0 < P, P <1)) ~= sum(and(0<S_prime, S_prime<1)))
    [m,n] = size(P);
    t = zeros(n,1);
    t(single_probe3(P,m,n)) = 1;
    r = zeros(m,1); r(find(S_prime * t)) = 1;
    P = update_probability(t, P, r);
    indices = NaN(1,n);
    
    
    for i = 1:n
        for j = 1:n_s
            if(sum(abs(P(:,i)-S(:,j))< 0.0000000000001) == m)
                indices(i) = i;
                P_prime(:,j) = P(:,i);
            end
        end
    end
    P_prime;
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr1 = n_itr1+1;
end
toc

%% Bundle Probing %%

S = diag(ones(96,1));    % Diagonal structure
S(:,96) = 1; S(96,:) = 1;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 20; batch_size = round(n_s/K);
P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
P_prime2 = zeros(m_s, n_s);
%P = ones(m_s, n_s)/20;
n_itr2 = 0;
tic
while( sum(and(0 < P, P <1)) ~= sum(and(0<S_prime, S_prime<1)))
    
    P_prime = P;
    [m_p, n_p] = size(P);
    probe_mat = zeros(n_p, K);r_mat = zeros(m_p, K);
    indices = NaN(1,n_p); ind = 1:n_p;
    % Finding the bundle of probes
    for i = 1:K
        [m,n] = size(P_prime);
        tau = single_probe3(P_prime,m,n);
        probe_mat(tau(1:min(size(tau,2),batch_size)),i) = 1;
        %ind = setdiff(ind, ind(tau(1:min(size(tau,2),batch_size))));
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
    % P_prime2
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr2= n_itr2+1;
end
toc

%% My METHOD %%
S = diag(ones(60,1));    % Diagonal structure
S(:,6) = 1; S(6,:) = 1;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 10; batch_size = round(n_s/K);
P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
P_prime2 = zeros(m_s, n_s);

n_itr3 = 0; probe_matrices = cell(5,1);
while( sum(and(0 < P, P <1)) ~= sum(and(0<S_prime, S_prime<1)))
    
    [m_p, n_p] = size(P);
    probe_mat = cluster_probe(P, n_p, min(K, n_p)); 
    r_mat = zeros(m_p, K);
    indices = NaN(1,n_p); ind = 1:n_p;
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
    % P_prime2
    
    % Remove the fully identified columns from the probability matrix
    % before the next iteration
    col_indices = (1:n_p).*(isnan(indices));
    P = P(:,col_indices(col_indices>0));
    % Remove the corresponding columns from the sparsity matrix
    S_prime = S_prime(:,col_indices(col_indices>0));
    
    n_itr3= n_itr3+1;
    probe_matrices{n_itr3}= probe_mat;
end

    
    
 
 
 