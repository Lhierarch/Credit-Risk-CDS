%% Author: Xinyu Zhang
%% CID: 01787342
%% Firm values in structured credit model
%
% Jump Diffusion Option Valuation in Discrete Time 

% JDEuropeanInstructions(40,30,0.08,1,sqrt(0.05),0,5,sqrt(0.05),-0.5*0.05^2,200,2)
% We will simulate a european put option in this case, here are the input
% parameters used:


% S0        =   Current Price.
% X         =   Strike Price.
% r         =   Risk-Free rate.
% sigma     =   Volatility of the underlying (assume lognormal process).
% D         =   Dividend Yield.
% T         =   Time to Maturity.

% lambda    =   Intensity of jump process.
% v         =   S.D. of jump size
% mu        =   Mean of jump size author suggests this as -0.5*v^2 because
% of how it modifies K (See next cell)

% flag1     =   1 for call 2 for put.
% n         =   Number of time steps.




%%

results_i=zeros(20,4);
V_i=[80,90,125,180];
T_i=0.5:0.5:10;

for i=1:20
    for j=1:4
        %%
        sigma=0.23; lambda=0.1; a=log(0.9); k=exp(a)-1; fi=0.4; rf=0.028; D=95;
        lambda2=lambda*(1+k);
        V=V_i(j); T=T_i(i); Call=0;
        for n=0:5
            v2=sigma^2+n*(fi^2)/T;
            rn=rf-lambda*k+n*a/T;
            % CallBS=Q1iii_fx(S0,X,r,T,sigma,D,lambda,v,mu,n,flag1)
            CallBS=Q1iii_fx(V,T,D,v2,rn);
            Call=Call+CallBS*exp(-lambda2*T)*((lambda2*T)^n)/factorial(n);
        end
        B=V-Call;
        %B=D*exp(-rf*T)-PutBS;
        y=(log(D/B))/T;
        s=y-rf;
        results_i(i,j)=s;
        
        
    end
end
writematrix(results_i,'results_i.xlsx','Sheet',1);

xlabel('Maturity t(years)')
ylabel('Credit Spread (b.p.)')
x=T_i;
y1_i=results_i(:,1);
y2_i=results_i(:,2);
y3_i=results_i(:,3);
y4_i=results_i(:,4);

y1_a=results_a(:,1);
y2_a=results_a(:,2);
y3_a=results_a(:,3);
y4_a=results_a(:,4);

plot(x,y1_i,x,y2_i,x,y3_i,x,y4_i,x,y1_a,x,y2_a,x,y3_a,x,y4_a)
