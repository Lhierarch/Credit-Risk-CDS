function debtfv = Q1i_fx(K)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rf=0.031;T=20;dt=1/12;n=T/dt;volatility1=0.15; volatility2=0.25; 
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

debtfv=abs(PD1-PD2);
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
end

