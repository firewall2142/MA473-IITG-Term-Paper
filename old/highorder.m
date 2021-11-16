function [Ss, V] = highorder(n,m,Smin,Smax,sigma,X,r,T,chi)
    % n -> number of space steps
    % m -> number of time steps
    % X -> strike price
    % returns Ss -> stocks points, V -> value at tau
    
    %%%
    % Tests
    %{
    (n,m,Smin,Smax,sigma,X,r,T,chi)
    highorder(10,10,0,45,0.3,15,0.02,0.5,12)
    %}
    %%%
    
    c1 = asinh(chi*(Smin-X)); c2 = asinh(chi*(Smax-X));
    psi = @(S) (asinh(chi*(S-X))-c1)/(c2-c1); % y = psi(S)
    phi = @(y) (1/chi)*sinh(c2*y + c1*(1-y)) + X; % S = phi(y)
    J = @(y) (1/chi)*((c2-c1)*cosh(c2*y + c1*(1-y))); %phi'(y)
    H = @(y) (1/chi)*((c2-c1)^2)*sinh(c2*y + c1*(1-y)); %phi''(y)
    
    ymin = psi(Smin); ymax = psi(Smax);
    h = (ymax-ymin)/n;
    ys = ((1:n)*(ymax-ymin)/n + ymin)';
    
    %%% 
    % Initialize D, D2 operators
    
    G = zeros(n);
    G(1,1)=-17/(6*h); G(1,2)=3/(2*h); G(1,3)=3/(2*h); G(1,4)=-1/(6*h);
    G(n,n)=17/(6*h);G(n,n-1)=-3/(2*h);G(n,n-2)=-3/(2*h);G(n,n-3)=1/(6*h);
    for i=2:n-1
        G(i,i-1) = -3/(4*h);
        G(i,i+1) =  3/(4*h);
    end
    
    F = eye(n);
    F(1,2) = 3; F(n,n-1) = 3;
    for i=2:n-1
        F(i,i-1) = 1/4;
        F(i,i+1) = 1/4;
    end
    
    W = zeros(n);
    W(1,1)=145/(12*h*h); W(1,2)=-76/(3*h*h);
    W(1,3)=29/(2*h*h); W(1,4)=-4/(3*h*h);W(1,5)=1/(12*h*h);
    
    W(n,n)=145/(12*h*h);W(n,n-1)=-76/(3*h*h);W(n,n-2)=29/(2*h*h);
    W(n,n-3)=-4/(3*h*h); W(n,n-4)=1/(12*h*h);
    
    for i=2:n-1
        W(i,i-1) = 6/(5*h*h);
        W(i,i) = -12/(5*h*h);
        W(i,i+1) = 6/(5*h*h);
    end
    
    U = eye(n);
    U(1,2) = 10; U(n,n-1) = 10;
    for i=2:n-1
        U(i,i-1) = 1/10;
        U(i,i+1) = 1/10;
    end
    
    D = F\G; % faster way of calculating inv(F) * G
    D2 = U\W;

    %%% 
    
    L = 0.5*(sigma^2)*diag((phi(ys)./J(ys)).^2)*D2 + ...
        diag(r*phi(ys)./J(ys) - ...
            0.5*(sigma^2)*(((phi(ys).^2)./(J(ys).^3)).*H(ys)) ...
        )*D - r*eye(n);
    
    del = T/m;
    ts = (0:m)*del;
    
    V = zeros(n,m+1);
    V(:,1) = max(phi(ys)-X,0);
    V(:,2) = (eye(n)-del*L/2)\((eye(n)+del*L/2)*V(:,1));
    V(:,3) = (eye(n)-del*L/2)\((eye(n)+del*L/2)*V(:,2));
    V(:,4) = (11*eye(n)/6 - del*L)\(3*V(:,3) - 1.5*V(:,2) + (1/3)*V(:,1));
    
    for j=5:m+1
        V(:,j) = (25*eye(n)/12 - del*L) \ ...
            (4*V(:,j-1)-3*V(:,j-2)+(4/3)*V(:,j-3)-0.25*V(:,j-4));
    end
    
    Ss = phi(ys);
end

