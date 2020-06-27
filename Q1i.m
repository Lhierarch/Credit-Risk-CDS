%% Author: Xinyu Zhang
%% CID: 01787342
%% Firm values in structured credit model
%
% This example shows how to build credit models.
%
% One of the fundamental tasks in credit risk management is to 
% compute the value of the firm (V ) for a 20-year period using the Merton
% model. 
% Use the following inputs:
% Risk-free rate = 3.10% 
% Face value of debt at maturity equals K = 80
% T = 20 and frequency is monthly (i.e. dt = 1/12)
% Asset volatility 15% for the first 10 years and 25% for the last 10 years
% V(0)=100.

rf=0.031;
K=80;
T=20;
dt=1/12;
n=T/dt;
volatility1=0.15; 
volatility2=0.25; 
firmvalue=zeros(250000,n);
firmvalue(:,1)=100;
ramdon=randn(250000,n);

for j=2:n/2
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility1*sqrt(dt);
end

for j=(n/2+1):n
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility2*sqrt(dt);
end

%%
%%(a)
% Simulate the value of the firm (using the above assumptions) 250,000 times. 
% Plot 5 randomly selected paths
x=1:n;
y1=firmvalue(1,:);
y2=firmvalue(2,:);
y3=firmvalue(3,:);
y4=firmvalue(4,:);
y5=firmvalue(5,:);
xlabel('Time t(months)')
ylabel('Value of firms V(t)')
title(' Simulated Asset Values ')
plot(x,y1,x,y2,x,y3,x,y4,x,y5)

%%
%%(b)
residualfv=firmvalue(:,n)-K;
residualfv=round(residualfv,0);
tbl=tabulate(residualfv);
plot(tbl(1:500,1),tbl(1:500,3)/100);

%%
%%c
sum1=0;
for i=1:250000
    if firmvalue(i,n)<K
        sum1=sum1+1;
    end
end
PD1=sum1/250000;
% PD1
% PD =

%    0.3266
mertonmodel

[PD,DD,A,Sa] = mertonmodel(Equity,EquityVol,Liability,Rate,'Drift',Drift)


%%
%%d
C=zeros(1,n);
for i=1:n
    C(i)=K/exp(rf*(n-i)*dt);
end

sum2=0;
for i =1:250000
    for j=1:n
        if firmvalue(i,j)<C(j)
            sum2=sum2+1;
            break;
        end
    end
end
PD2=sum2/250000;
% PD2

%%
%e
%the default probabilities assuming model 1 and model 2 significantly different? 
% I can run the codes before several times
rf=0.031;K=80;T=20;dt=1/12;n=T/dt;volatility1=0.15; volatility2=0.25; 
firmvalue=zeros(250000,n);firmvalue(:,1)=100;ramdon=randn(250000,n);
for j=2:n/2
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility1*sqrt(dt);
end
for j=(n/2+1):n
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility2*sqrt(dt);
end
sum1=0;
for i=1:250000
    if firmvalue(i,n)<K
        sum1=sum1+1;
    end
end
PD1=sum1/250000;
C=zeros(1,n);
for i=1:n
    C(i)=K/exp(rf*(n-i)*dt);
end
sum2=0;
for i =1:250000
    for j=1:n
        if firmvalue(i,j)<C(j)
            sum2=sum2+1;
            break;
        end
    end
end
PD2=sum2/250000;
% later this part being put into the function part for Q1i
% to calculate the difference between PD1 and PD2
% PD1=0.3259  PD2=0.4995
% PD1=0.3246  PD2=0.4979
% PD1=0.3248  PD2=0.4984
% PD1=0.3263  PD2=0.4989
% PD1=0.3245  PD2=0.4973
% 

%%
%f
% x = fminsearch(fun,x0,options) minimizes with the optimization options 
% specified in the structure options. Use optimset to set these options

options = optimset('PlotFcns','optimplotfval','TolFun',1e-7);
fun=@Q1i_fx; % my own function
x0=[100];
[x,fval] = fminsearch(fun,x0,options);

%x=30;fval=0.0373

%%
%%g
rf=0.031;K=80;T=20;
dt=1/1;  % this time dt change to 1/1
n=T/dt;volatility1=0.15; volatility2=0.25; 
firmvalue=zeros(250000,n);firmvalue(:,1)=100;ramdon=randn(250000,n);
for j=2:n/2
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility1*sqrt(dt);
end
for j=(n/2+1):n
    firmvalue(:,j)=firmvalue(:,j-1)+firmvalue(:,j-1).*rf*dt+firmvalue(:,j-1).*ramdon(:,j).*volatility2*sqrt(dt);
end
sum1=0;
for i=1:250000
    if firmvalue(i,n)<K
        sum1=sum1+1;
    end
end
PD1=sum1/250000;
C=zeros(1,n);
for i=1:n
    C(i)=K/exp(rf*(n-i)*dt);
end
sum2=0;
for i =1:250000
    for j=1:n
        if firmvalue(i,j)<C(j)
            sum2=sum2+1;
            break;
        end
    end
end
PD2=sum2/250000;
% PD1       PD2
% 0.3346    0.4557

pd=makedist("normal","mu",0,"sigma",1);
y=cdf(pd,-0.3166);
%y=0.3758

