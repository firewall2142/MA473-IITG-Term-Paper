function [V] = BDF4(n,m,T,V0,A,b)
%BDF4(n,m,T,V0,A,b) calculates option prices V : (n,m+1) matrix
%   V0 is the option price at Maturity
%   T is the time to maturity
%   n -> space dim, m -> time dim
    V = zeros(n,m+1);
    del = T/m;
    V0 = reshape(V0,n,1);
    V(:,1) = V0;
    V(:,2) = (eye(n)-del*A/2)\((eye(n)+del*A/2)*V(:,1)+del*b);
    V(:,3) = (eye(n)-del*A/2)\((eye(n)+del*A/2)*V(:,2)+del*b);
    V(:,4) = (11*eye(n)/6 - del*A)\(3*V(:,3) - 1.5*V(:,2) + ...
        (1/3)*V(:,1) + del*b);
    
    for j=5:m+1
        V(:,j) = (25*eye(n)/12 - del*A) \ ...
            (4*V(:,j-1)-3*V(:,j-2)+(4/3)*V(:,j-3)-0.25*V(:,j-4)+del*b);
    end
end

