%%%% Generalize the Algorithm %%%%

%    Define the sparsity matrix

S = diag(ones(4,1));    % Diagonal structure
S(:,4) = 1; S(4,:) = 1; % Arrow type structure
 
%%% Sparsity matrix from the suitesparse.com website
%   Weblink: http://yifanhu.net/GALLERY/GRAPHS/search.html

load cage4.mat
sparsemat = full(Problem.A); 
S = abs(sign(sparsemat)); % The sparse matrix

[m_s, n_s] = size(S);     % Size of the sparse matrix
S_prime = S;
% Define the probability matrix
P = ones(m_s,n_s)/min(m_s,n_s); % We begin with P matrix with identical elements 1/min(m, n)
P_prime = zeros(m_s,n_s); 
n_itr =0;                 % Counter for the number of iterations


 while( sum(and(0 < P, P <1)) ~= sum(and(0<S_prime, S_prime<1)))
     [m,n] = size(P); 
     t = zeros(n,1);
     t(single_probe3(P,m,n)) = 1;
     r = zeros(m,1); r(find(S_prime * t)) = 1;
     P = update_probability(t, P, r);
     indices = NaN(1,n);
     
     
     for i = 1:n
         for j = 1:n_s
             if(sum(abs(P(:,i)-S(:,j))< 0.0000001) == m)
             indices(i) = i;
             P_prime(:,j) = P(:,i);
             end
         end
     end
     
     
     % Remove the fully identified columns from the probability matrix
     % before the next iteration
     col_indices = (1:n).*(isnan(indices));
     P = P(:,col_indices(col_indices>0));
     % Remove the corresponding columns from the sparsity matrix
     S_prime = S_prime(:,col_indices(col_indices>0));
     
     n_itr = n_itr+1;
 end
 P_prime
 n_itr
 
  %% Bundle Probing %%
 

S = diag(ones(160,1));    % Diagonal structure
S(:,160) = 1; S(160,:) = 1;
S_prime = S;
[m_s, n_s] = size(S);     % Size of the sparse matrix
K = 32; batch_size = round(n_s/K);
P = diag(sum(S,2)/n_s)*ones(m_s, n_s); 
P_prime2 = zeros(m_s, n_s);
%P = ones(m_s, n_s)/20;
n_itr2 = 0;

while( sum(and(0 < P, P <1)) ~= sum(and(0<S_prime, S_prime<1)))
   
P_prime = P;
[m_p, n_p] = size(P);
probe_mat = zeros(n_p, K);r_mat = zeros(m_p, K);
indices = NaN(1,n_p);
% Finding the bundle of probes
for i = 1:K
    [m,n] = size(P_prime);  
    tau = single_probe3(P_prime,m,n);
    probe_mat(tau(1:min(size(tau,2),batch_size)),i) = 1;
    P_prime(:,tau(1:min(size(tau,2),batch_size)))= ones(m,min(size(tau,2),batch_size)) ;
end
% The corresponding result matrix
r_mat(find(S_prime * probe_mat)) = 1;

    P = update_probability(probe_mat, P, r_mat);

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
 
