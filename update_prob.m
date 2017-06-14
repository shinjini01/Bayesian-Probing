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
n_itr =0;                 % Counter for the number of interations


 while( sum(sum(0<P<1)) ~=sum(sum(0<S_prime<1)))
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
     P_prime
     
     % Remove the fully identified columns from the probability matrix
     % before the next iteration
     col_indices = (1:n).*(isnan(indices));
     P = P(:,col_indices(col_indices>0));
     % Remove the corresponding columns from the sparsity matrix
     S_prime = S_prime(:,col_indices(col_indices>0));
     
     n_itr = n_itr+1;
 end
 
