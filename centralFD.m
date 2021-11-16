function [P,q] = centralFD(n,Smin,Smax,sigma,X,r)
%CENTRALFD(n,Smin,Smax,sigma,X,r)
%   returns P,q s.t. dV/dt = PV + q

    h = (Smax-Smin)/(n+1);
    A = zeros(n);
    A(1,1:4) = [-10 18 -6 1];
    A(n,n-3:n) = [1 -6 18 -10];
    for i=2:n-1
        A(i,i-1:i+1) = [-8 0 8];
        if i-2>0
            A(i,i-2)=1;
        end
        if i+2<=n
            A(i,i+2)=-1;
        end
    end
    A = A/(12*h);
    Vn = Smax-X;
    b1 = zeros(n,1);
    b1(n-1:n) = [-Vn -3*Vn]/(12*h);
    
    B = zeros(n);
    B(1,1:5)=[-15 -4 14 -6 1];
    B(n,n-4:n)=[1 -6 14 -4 -15];
    for i=2:n-1
        B(i,i-1:i+1) = [16 -30 16];
        if i-2>0
            B(i,i-2)=-1;
        end
        if i+2<=n
            B(i,i+2)=-1;
        end
    end
    B = B/(12*h*h);
    b2 = zeros(n,1);
    b2(n-1:n) = [-Vn  10*Vn]/(12*h*h);
    
    S = linspace(Smin,Smax,n+2);
    S = S(2:n+1)';
    P = -((1/2)*sigma^2*(S.^2).*B + r*S.*A - r*eye(n));
    q = -((1/2)*sigma^2*(S.^2).*b2 + r*S.*b1);
end
