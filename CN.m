function [V] = CN(n,m,T,V0,A,b)
%CN(n,m,T,V0,A,b) calculates option prices V : (n,m+1) matrix
%   V0 is the option price at Maturity
%   T is the time to maturity
%   n -> space dim, m -> time dim
    V = zeros(n,m+1);
    del = T/m;
    V0 = reshape(V0,n,1);
    V(:,1) = V0;
    
    for j=2:m+1
        V(:,j) = (eye(n)-del*A/2)\((eye(n)+del*A/2)*V(:,j-1)+del*b);
    end
end

