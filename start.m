clear


%% Compact + BDF4

Smin=0; Smax=45;
sigma=0.3; X=15; r=0.02;
T=0.5;
Ns = [10 20 40];

errs = [];
for n=Ns
    m = n;
    [A,b]=compactFD(n,Smin,Smax,sigma,X,r);
    Ss = linspace(Smin,Smax,n+2);
    Ss = Ss(2:n+1);
    V0 = max(Ss-X, 0);
    V = BDF4(n,m,T,V0,A,b);
    Vt = eurocall(Ss,sigma,X,r,T)';
    errs = [errs mean(abs(V(:,end)-Vt),'all')];
end
ords = errs*0;
for i=2:length(errs)
    ords(i) = -log(errs(i-1)/errs(i))/log(Ns(i-1)/Ns(i));
end

fprintf("Grid Size\tError\tOrder\n");
for i=1:length(Ns)
    n = Ns(i);
    err = errs(i);
    ord = ords(i);
    fprintf("%dx%d\t%f\t%f\n",n,n,err,ord);
end


figure;
hold on
plot(Ss,V(:,end),'-*');
plot(Ss,eurocall(Ss,sigma,X,r,T));
legend('Calculated','Theoretical');
xlabel('S');
ylabel('V');
saveas(gcf,'output/comp_bdf4.png');

figure;
plot(Ss,abs(V(:,end)-eurocall(Ss,sigma,X,r,T)'));
xlabel('S');
ylabel('Error');
saveas(gcf,'output/comp_bdf4_e.png');

%% Compact + CN
Smin=0; Smax=45;
sigma=0.3; X=15; r=0.02;
T=0.5;
Ns = [10 20 40];

errs = [];
for n=Ns
    m = n;
    [A,b]=compactFD(n,Smin,Smax,sigma,X,r);
    Ss = linspace(Smin,Smax,n+2);
    Ss = Ss(2:n+1);
    V0 = max(Ss-X, 0);
    V = CN(n,m,T,V0,A,b);
    Vt = eurocall(Ss,sigma,X,r,T)';
    errs = [errs mean(abs(V(:,end)-Vt),'all')];
end
ords = errs*0;
for i=2:length(errs)
    ords(i) = -log(errs(i-1)/errs(i))/log(Ns(i-1)/Ns(i));
end

fprintf("Grid Size\tError\tOrder\n");
for i=1:length(Ns)
    n = Ns(i);
    err = errs(i);
    ord = ords(i);
    fprintf("%dx%d\t%f\t%f\n",n,n,err,ord);
end


figure;
hold on
plot(Ss,V(:,end),'-*');
plot(Ss,eurocall(Ss,sigma,X,r,T));
legend('Calculated','Theoretical');
xlabel('S');
ylabel('V');
saveas(gcf,'output/comp_cn.png');

figure;
plot(Ss,abs(V(:,end)-eurocall(Ss,sigma,X,r,T)'));
xlabel('S');
ylabel('Error');
saveas(gcf,'output/comp_cn_e.png');


%% CompactFD + CN + GridRefine
Smin=0; Smax=45;
sigma=0.3; X=15; r=0.02;
T=0.5;
Ns = [10 20 40];
xi = 12;

errs = [];
for n=Ns
    m = n;
    [A,b,Ss]=grCompFD(n,Smin,Smax,sigma,X,r,xi);
    V0 = max(Ss-X, 0);
    V = CN(n,m,T,V0,A,b);
    Vt = eurocall(Ss,sigma,X,r,T);
    errs = [errs mean(abs(V(:,end)-Vt),'all')];
end
ords = errs*0;
for i=2:length(errs)
    ords(i) = -log(errs(i-1)/errs(i))/log(Ns(i-1)/Ns(i));
end

fprintf("Grid Size\tError\tOrder\n");
for i=1:length(Ns)
    n = Ns(i);
    err = errs(i);
    ord = ords(i);
    fprintf("%dx%d\t%f\t%f\n",n,n,err,ord);
end


figure;
hold on
plot(Ss,V(:,end),'-*');
plot(Ss,eurocall(Ss,sigma,X,r,T));
legend('Calculated','Theoretical');
xlabel('S');
ylabel('V');
saveas(gcf,'output/comp_cn_gr.png');

figure;
plot(Ss,abs(V(:,end)-eurocall(Ss,sigma,X,r,T)));
xlabel('S');
ylabel('Error');
saveas(gcf,'output/comp_cn_gr_e.png');


%% CompactFD + BDF4 + GridRefine
Smin=0; Smax=45;
sigma=0.3; X=15; r=0.02;
T=0.5;
Ns = [10 20 40];
xi = 12;

errs = [];
for n=Ns
    m = n;
    [A,b,Ss]=grCompFD(n,Smin,Smax,sigma,X,r,xi);
    V0 = max(Ss-X, 0);
    V = BDF4(n,m,T,V0,A,b);
    Vt = eurocall(Ss,sigma,X,r,T);
    errs = [errs mean(abs(V(:,end)-Vt),'all')];
end
ords = errs*0;
for i=2:length(errs)
    ords(i) = -log(errs(i-1)/errs(i))/log(Ns(i-1)/Ns(i));
end

fprintf("Grid Size\tError\tOrder\n");
for i=1:length(Ns)
    n = Ns(i);
    err = errs(i);
    ord = ords(i);
    fprintf("%dx%d\t%f\t%f\n",n,n,err,ord);
end


figure;
hold on
plot(Ss,V(:,end),'-*');
plot(Ss,eurocall(Ss,sigma,X,r,T));
legend('Calculated','Theoretical');
xlabel('S');
ylabel('V');
saveas(gcf,'output/comp_bdf4_gr.png');

figure;
plot(Ss,abs(V(:,end)-eurocall(Ss,sigma,X,r,T)));
xlabel('S');
ylabel('Error');
saveas(gcf,'output/comp_bdf4_gr_e.png');


