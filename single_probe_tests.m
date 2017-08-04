
 S = diag(ones(4,1));    % Diagonal structure
 S = speye(1000);
 S(:,10) = 1; S(10,:) = 1;
 [m_s, n_s] = size(S); 

load ibm32.mat
sparsemat = full(Problem.A);
S = abs(sign(sparsemat));
n_itr2 = 0;
n_itr2 = zeros(3,1);
for l = 1:3
    
    [m_s, n_s] = size(S);     % Size of the sparse matrix
  
    if(l ==1)
        P = single(ones(m_s, n_s)/n_s);
    elseif (l==2)
        P = diag(sum(S,2)/n_s)*ones(m_s, n_s);
    elseif (l==3)
        P = (sum(sum(S))/(m_s*n_s))*ones(m_s, n_s);
        %P = (3/n_s)*ones(m_s, n_s);
    end
   
    n_itr2(l) = find_probe(P,S);
end

n_itr2 = zeros(3,1);
for l = 1:3
        % Size of the sparse matrix
  
    if(l ==1)
        n_itr = 0;
        for i = 1:100
            P = rand(n_s)*(1 - (2/n_s)) + (1/n_s);
            n_itr = n_itr + find_probe(P,S);
        end
        n_itr2(l) = n_itr/100;
        elseif (l==2)
            n_itr = 0;
            for i = 1:100
                mnzi = rand(m_s, n_s)*(2*log(n_s) -1) + 1;nnzj = rand(m_s, n_s)*(2*log(n_s) -1) + 1;
                P = min(mnzi,nnzj)/n_s;
                n_itr = n_itr + find_probe(P,S);
            end
            n_itr2(l) = n_itr/100;
            elseif (l==3)
                n_itr =0;
                for i = 1:100
                    a = random('beta', (1/n_s), (2-1/n_s), m_s, n_s, 1000);
                    b = random('beta', (1+1/n_s), (1-1/n_s), m_s, n_s, 1000);
                    P = (mean(a, 3)+mean(b,3))/2;
                    n_itr = n_itr + find_probe(P,S);
            end
            n_itr2(l) = n_itr/100;
    end
end
S = diag(ones(100,1)); 
S = speye(10000);
[m_s, n_s] = size(S); n_itr2 = 0; P = ones(m_s, n_s)/n_s;

