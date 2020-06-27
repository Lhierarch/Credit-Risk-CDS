function tcopula=Q3_tcopula(a)
rho=a;
N=250000;
%lambda=0.15;
lambda=[ones(N,10)*0.05,ones(N,10)*0.15,ones(N,10)*0.25];

sigma=ones(30,30).*rho;
for i=1:30
    sigma(i,i)=1;
end

mu=zeros(1,30);

rng default;   % For reproducibility% For reproducibility
R=mvtrnd(sigma,4,250000);
figure;
plot(R(:,1),R(:,2),'+')
%%

U=tcdf(R,4);

tau=-(1./lambda).*log(1-U);

xbar = mean(R);  % Sample mean
s = std(R);      % Sample standard deviation
t = (xbar - mu)/(s/sqrt(N))
p = 1-tcdf(t,N-1) % Probability of larger t-statistic

%%
r1=sum(sum(tau<=1.5,2)>=8)/N;
r2=sum(sum(tau<=1.5,2)>=10)/N;
r3=sum(sum(tau<=1.5,2)>=12)/N;
r4=sum(sum(tau<=2,2)==0)/N;
r5=sum(sum(tau<=3,2)==0)/N;
r6=sum(sum(tau<=3.5,2)==0)/N;
tcopula=[r1,r2,r3,r4,r5,r6];
end

