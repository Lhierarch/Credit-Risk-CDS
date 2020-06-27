function Call=Q1iii_fx(V,T,D,sigma2,rf)
%V,T,D,sigma2,rf

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

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%K=exp(mu+0.5*v^2)-1;                % Expectation(Y-1). (relative expected jump size)
%EY=exp(mu+1/2*v^2);                 % Expectation of y (absolute expected jump size)

%h=T/n;                              % Time Step (dt)
%alpha=r-D-1/2*sigma^2-lambda*K;     % Drift

%%a=(exp(r*h-D*h)-lambda*h*EY)/(1-lambda*h);   % (drift minus expected jump size during period h)/(probability of jump not occurring during h) NEEDS TO BE EY
%b=exp(alpha*h-sigma*sqrt(h));                % d
%c=exp(alpha*h+sigma*sqrt(h));                % u

%q=(a-b)/(c-b)


pd=makedist("normal","mu",0,"sigma",1);
    sigma=sqrt(sigma2);
    d1=(log(V/D)+(rf+0.5*sigma^2)*T)/(sigma*sqrt(T));
    d2=d1-sigma*sqrt(T);
    Nd1=cdf(pd,d1);
    Nd2=cdf(pd,d2);
    Call=V*Nd1-D*Nd2*exp(-rf*T);
    
end

