n=50;m=50;
Smin=0; Smax=45;
sigma=0.3; X=15; r=0.02;
T=0.5;
chi=12;

[Ss,V1] = highorder(n,m,Smin,Smax,sigma,X,r,T,chi);
V2 = eurocall(Ss,sigma,X,r,T);

tiledlayout(2,2);

nexttile
plot(Ss, V1(:,end),'-*');
title('High order scheme');

nexttile
plot(Ss, V2,'-*');
title('True solution');

nexttile
plot(Ss, abs(V1(:,end)-V2)/V2);
title('Relative error');

nexttile
plot(Ss, abs(V1(:,end)-V2));
title('Absolute error');