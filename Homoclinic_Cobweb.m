%% Homoclinic Cobweb
% This script calculates the unstable manifold for a cobweb model with 
% heterogenous expectations as described by Hommes (2013, Ch. 5) which 
% for specific parameter values demonstrates a homoclinic bifurcation
% (Written by Kian Abbas Nejad)

%% Housekeeping
clearvars
clc
close all

tmax=40000;
trans = 5000; %number of transient iterations
p = zeros(1,tmax); % price
u1 = zeros(1,tmax); %Utility of type 1
u2 = zeros(1,tmax); % utility of type 2
n1 = zeros(1,tmax); %fraction of type 1
n2 = zeros(1,tmax); % fraction of type 2
h1 = zeros(1,tmax); % expectation of type 2
h2 = zeros(1,tmax); % expectation of type 2

% Parameters
beta = 3;
delta = 0.5; % strategy updating parameter
beta1 = -0.85; %forecasting rule for contrarian
c1 = 1; 
c2 = 0;
a = 0;
d = 0.5;
s = 1.35;

%Initial Conditions -- need 3 boundary conditions
p(1)=2;
n1(1)=0.5;
n2(1)=0.5;
pstar = a/(d+s); % RE price

% Formulae 
fu1= @(T,P,H1) s*P(T)*H1(T)-(s/2)*(H1(T)).^2-c1; 
fu2= @(T,P,H2) s*P(T)*H2(T)-(s/2)*(H2(T)).^2-c2; 
fn1 = @(T,U1,U2,Beta,N1) (1-delta)*exp(Beta*U1(T))/(exp(Beta*U1(T))+exp(Beta*U2(T)))+delta* N1(T-1);
fn2 = @(T,U1,U2,Beta,N2) (1-delta)*exp(Beta*U2(T))/(exp(Beta*U1(T))+exp(Beta*U2(T)))+delta* N2(T-1);
fp = @(T,N1,N2,H1,H2) (a-N1(T-1)*s*H1(T)-N2(T-1)*s*H2(T))./d;

% specific types
fh1 = @(P,T) pstar + beta1*(P(T-1)-pstar);
fh2 = @(P,T) P(T-1); 


% ODE System Starts at period 2
for i=2:tmax
    h1(i)= fh1(p,i);
    h2(i)= fh2(p,i);
    p(i) = fp(i,n1,n2,h1,h2);
    u1(i) = fu1(i,p,h1);
    u2(i) = fu2(i,p,h2);
    n1(i) = fn1(i,u1,u2,beta,n1);
    n2(i) = fn2(i,u1,u2,beta,n2);
end

%Attractor
fig1=figure(1);
x = p - pstar;
scatter(x(trans:end),n1(trans:end),0.7,'k')
fig1.Name='Attractor';
title('strange attractor for $\beta=3$','FontSize',14,'interpreter','latex');
ylabel('$n_{1,t}$','FontSize',16,'interpreter','latex');
xlabel('$x_t$','FontSize',16,'interpreter','latex');


%% Unstable Manifold

tmax = 10; % number of iterations
tmax = tmax+1;
len = 10000; %length of initial vector
dx = 0.001; % perturbation size

 p = zeros(1,tmax); 
 u1 = zeros(1,tmax); 
 u2 = zeros(1,tmax); 
 n1 = zeros(1,tmax); 
 n2 = zeros(1,tmax); 
 h1 = zeros(1,tmax); 
 h2 = zeros(1,tmax); 

% steady state
n1ss = (tanh(-beta*c1/2)+1)/2;

%Initial vector
pvec = zeros(len,tmax);
n1vec = zeros(len,tmax);
pvec(:,1) = linspace(pstar+dx,pstar-dx,len);
n1vec(:,1)= linspace(n1ss,n1ss,len);

% ODE System 
for k=1:len
    p(1)=pvec(k,1);
    n1(1)=n1vec(k,1);
    n2(1)=1-n2(1);
    for i=2:tmax
        h1(i)= fh1(p,i);
        h2(i)= fh2(p,i);
        p(i) = fp(i,n1,n2,h1,h2);
        u1(i) = fu1(i,p,h1);
        u2(i) = fu2(i,p,h2);
        n1(i) = fn1(i,u1,u2,beta,n1);
        n2(i) = fn2(i,u1,u2,beta,n2);
    end
    pvec(k,:)=p;
    n1vec(k,:)=n1;
end

fig2=figure(2);
hold on
fig2.Name='Manifold';
title('$\beta=3$: transversal homoclinic self-intersections','FontSize',14,'interpreter','latex');
ylabel('$n_{1,t}$','FontSize',16,'interpreter','latex');
xlabel('$x_t$','FontSize',16,'interpreter','latex');
for i=1:tmax
    scatter(pvec(:,i),n1vec(:,i),1,'k')
end





