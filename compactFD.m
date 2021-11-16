function [L,b] = compactFD(n,Smin,Smax,sigma,X,r)
%CENTRALFD(n,Smin,Smax,sigma,X,r)
%   returns P,q s.t. dV/dt = PV + q
    
    h = (Smax-Smin)/(n+1);
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
    
    S = linspace(Smin,Smax,n+2)';
    S = S(2:n+1);
    
    L = 0.5*(sigma^2)*(S.^2).*D2 + r*S.*D - r*eye(n);
    b = zeros(n,1);
end
