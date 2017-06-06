P = ones(8,8)/8 ;    % Initialize the probability matrix
Q = 1 - P ;          % 1 - Probability matrix, required to calculate w values
% To initialize the matrix with each column as a probe vector
a1 = [ones(4,1); zeros(4,1)];
a2 = repmat([1;1;0;0],2,1);
a3 = repmat([1;0],4,1);
t = [a1, 1-a1, a2, 1 - a2, a3, 1 - a3]; % t is the required matrix of probes
a0 = repmat([1;0],3,1); 
a4 = [a0; 1;1]; 
a5 = [1-a0; 1;1];
r.mat = [a1, 1 - a1, a2, 1 - a2, a4, a5]; % r.mat is the result matrix



for k = 1:size(t,2)
    P = update_probability(t(:,k), P, r.mat(:,k))
    end
    
