function GaussianCopula=Q3_Gcopula(a)

rho=a;
N=250000;
%lambda=0.15;
%
lambda=[ones(N,10)*0.05,ones(N,10)*0.15,ones(N,10)*0.25];

% R = mvnrnd(mu,sigma,n) returns a matrix R of n random vectors chosen 
% from the same multivariate normal distribution, with mean vector mu 
% and covariance matrix sigma. 

sigma=ones(30,30).*rho;
% SIGMA = [1.00 , 0.80 ; 0.80 , 1.00 ] ; 
for i=1:30
    sigma(i,i)=1;
end
mu=zeros(1,30);
rng('default')  % For reproducibility
R=mvnrnd(mu,sigma,250000);

plot(R(:,1),R(:,2),'+') %plot the ramdon numbers to see

U=normcdf(R);
tau=-(1./lambda).*log(1-U);

r1=sum(sum(tau<=1.5,2)>=8)/N;
r2=sum(sum(tau<=1.5,2)>=10)/N;
r3=sum(sum(tau<=1.5,2)>=12)/N;
r4=sum(sum(tau<=2,2)==0)/N;
r5=sum(sum(tau<=3,2)==0)/N;
r6=sum(sum(tau<=3.5,2)==0)/N;

GaussianCopula=[r1,r2,r3,r4,r5,r6];


end

