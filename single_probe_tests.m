%Code to find the number of iterations required for convergence
%for different choices of the sparse matrix and different initial P matrix
%as discussed in section 1.4 of the report

%Choice of the sparse matrix
S = diag(ones(100,1));     % Diagonal structure
S = diag(ones(100,1)); S(:,100) = 1; S(100,:) = 1; % Arrow head structure

load can96.mat;            % Banded 96x96 matrix
load ibm32.mat;            % HB 32x32 matrix
sparsemat = full(Problem.A);
S = abs(sign(sparsemat));


[m_s, n_s] = size(S);      % Size of the sparse matrix

% Determining # of iterations for the first two choices of P matrix
n_itr = zeros(2,1); %Number of iterations required for each of the 2 choices
for l = 1:2
    
    if(l ==1)
        P = ones(m_s, n_s)/n_s;
    elseif (l==2)
        P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
    end
    
    n_itr(l) = find_probe(P,S);
end
% Determining # of iterations for the last four choices of P matrix
n_itr2 = zeros(4,1); %Number of iterations required for each of the 4 choices
for l = 1:4
    
    if(l ==1)
        %P initialized with uniform random samples
        n_itr2 = 0;
        for i = 1:100
            P = rand(n_s)*(1 - (2/n_s)) + (1/n_s);
            n_itr2 = n_itr2 + find_probe(P,S);
        end
        n_itr2(l) = n_itr2/100; %Average of 100 trials
    elseif (l==2)
        %P  = min(z)/n, z being random sample from (1, 2logn)
        n_itr2 = 0;
        for i = 1:100
            mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
            P = min(mnzi,nnzj)/n_s;
            n_itr2 = n_itr2 + find_probe(P,S);
        end
        n_itr2(l) = n_itr2/100;
    elseif (l==3)
        %P initialized using Beta(1+1/n, 1-1/n) samples
        n_itr2 =0;
        for i = 1:100
            
            a = random('beta', (1+1/n_s), (1-1/n_s), m_s, n_s, 500);
            P = mean(a, 3);
            n_itr2 = n_itr2 + find_probe(P,S);
        end
        n_itr2(l) = n_itr2/100;
    elseif (l==4)
        %P initialized using average of Beta(1/n, 2-1/n) and Beta(1+1/n, 1-1/n)samples
        n_itr2 =0;
        for i = 1:100
            a = random('beta', (1/n_s), (2-1/n_s), m_s, n_s, 1000);
            b = random('beta', (1+1/n_s), (1-1/n_s), m_s, n_s, 1000);
            P = (mean(a, 3)+mean(b,3))/2;
            n_itr2 = n_itr2 + find_probe(P,S);
        end
        n_itr2(l) = n_itr2/100;
    end
end

