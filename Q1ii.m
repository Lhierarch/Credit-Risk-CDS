%% Author: Xinyu Zhang
%% CID: 01787342
%% Merton Model Credit Spreads:
%
% Study the credit spread as a function of the time to maturity(0 to 10 years) 
% for V (value of firm) equal to 75, 90, 125, and 180

% Use the following inputs:
% Risk-free rate = 2.80%, Asset volatility = 23%.
% Face value of debt equals to K = 95
% Asset volatility = 23%.


%% RESHAPE VARIABLES TO 3D TO FACILITATE MATRIXWISE BIVAR NEWTON-REPHSON
[n,m]=size(E_t);
nm=n*m;
E_t=reshape(E_t,1,1,nm);sig_E=reshape(sig_E,1,1,nm); 
K=reshape(K,1,1,nm); T=reshape(T,1,1,nm);

%% DEFINE JACOBIAN ANONYMOUS FUNCTION FOR [A_t & sig_A] SOLUTION
d1=@(A_t,sig_A,C)((1./(sig_A(:,:,C).*sqrt(T(:,:,C)-t))).*(log(A_t(:,:,C)./K(:,:,C)) + (r + 0.5.*sig_A(:,:,C).^2).*(T(:,:,C)-t)));
d2=@(A_t,sig_A,C)((1./(sig_A(:,:,C).*sqrt(T(:,:,C)-t))).*(log(A_t(:,:,C)./K(:,:,C)) + (r - 0.5.*sig_A(:,:,C).^2).*(T(:,:,C)-t)));

% System of Nonlinear Equation
fcnF1=@(A_t,sig_A,C)((A_t(:,:,C).*fcnN(d1(A_t,sig_A,C))-K(:,:,C).*exp(-r*(T(:,:,C)-t)).*fcnN(d2(A_t,sig_A,C)))-E_t(:,:,C));
fcnF2=@(A_t,sig_A,C)(A_t(:,:,C).*sig_A(:,:,C).*fcnN(d1(A_t,sig_A,C))-sig_E(:,:,C).*E_t(:,:,C));
fcnF=@(A_t,sig_A,C)([fcnF1(A_t,sig_A,C);fcnF2(A_t,sig_A,C)]);

% Define Partial Derivative Functions of Equation Set
fcnJ11=@(A_t,sig_A,C)(fcnN(d1(A_t,sig_A,C)));
fcnJ12=@(A_t,sig_A,C)(A_t(:,:,C).*sqrt((T(:,:,C)-t)).*fcnn(d1(A_t,sig_A,C)));
fcnJ21=@(A_t,sig_A,C)(fcnN(d1(A_t,sig_A,C)).*sig_A(:,:,C)+(fcnn(d1(A_t,sig_A,C))./sqrt((T(:,:,C)-t))));
fcnJ22=@(A_t,sig_A,C)(A_t(:,:,C).*fcnN(d1(A_t,sig_A,C)) + A_t(:,:,C).*sig_A(:,:,C).*fcnn(d1(A_t,sig_A,C)).*((-log(A_t(:,:,C)./(K(:,:,C).*exp(-r.*(T(:,:,C)-t))))./((sig_A(:,:,C).^2).*sqrt(T(:,:,C)-t)))+(0.5*sqrt(T(:,:,C)-t))));


% Define Jacobian Matrix
fcnJ=@(A_t,sig_A,C)([fcnJ11(A_t,sig_A,C),fcnJ12(A_t,sig_A,C);fcnJ21(A_t,sig_A,C),fcnJ22(A_t,sig_A,C)]);


%% SOLVE FOR ASSET VALUE & VOLATILITY [A_t & sig_A]
tolMat=1e-10;
k_max=20;
k=1;

% Initial Estimates for A_t & sig_A
A_t=E_t+K;
sig_A=sig_E.*E_t./(E_t+K);
A_t=reshape(A_t,1,1,nm); sig_A = reshape(sig_A,1,1,nm); C = true(1,1,nm);
x=[A_t;sig_A];
while any(C) && k<=k_max
    dx=mtimesx(inv3d(fcnJ(x(1,:,:),x(2,:,:),C)),fcnF(x(1,:,:),x(2,:,:),C));
    x(:,:,C) = x(:,:,C) - dx;
    C=any(abs(fcnF(x(1,:,:),x(2,:,:),true(1,1,nm)))>tolMat,1);
    k=k + 1;
end
A_t=x(1,:,:); A_t = reshape(A_t,n,m);
sig_A=x(2,:,:); sig_A = reshape(sig_A,n,m);
E_t=reshape(E_t,n,m); K = reshape(K,n,m);  T = reshape(T,n,m); 

%% SOLVE FOR FIRM DEBT VALUE [D_t]
D_t=A_t - E_t; % Balance Sheet

%% SOLVE FOR CREDIT SPREAD [s]
s=-(log(D_t./K)./(T-t))-r;

%% SOLVE FOR DEFAULT PROBABILITY [p]
% Redefine Black-Scholes Functions Since Variables have been reshaped
d1=@(A_t,sig_A,C)((1./(sig_A(:,:,C).*sqrt(T(:,:,C)-t))).*(log(A_t(:,:,C)./K(:,:,C)) + (r + 0.5.*sig_A(:,:,C).^2).*(T(:,:,C)-t)));
d2=@(A_t,sig_A,C)((1./(sig_A(:,:,C).*sqrt(T(:,:,C)-t))).*(log(A_t(:,:,C)./K(:,:,C)) + (r - 0.5.*sig_A(:,:,C).^2).*(T(:,:,C)-t)));
p=fcnN(-d2(A_t,sig_A,true)); % P(A_t < K) = N(-d_m)

%% EXPECTED RECOVERY ON DEFAULT [R]
R = exp(r.*(T-t)) .* (A_t./K) .* (fcnN(-d1(A_t,sig_A,true))./fcnN(-d2(A_t,sig_A,true)));

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%MethodII%%%%%%%%%%%%%%%%%%%%%%%%%
%
%

%%
%a

V_a=[80,90,125,180];
sigma_a=0.23;
T_a=0.5:0.5:10;
rf_a=0.028;
results_a=zeros(20,4); %the table output here


for i=1:20
    for j=1:4
        
      %  do the loop for matrix
        V=V_a(j);
        T=T_a(i);
        rf=rf_a;  %  rf=rfAll(j);
        D=95; %  face value of debt
        sigma=sigma_a; %volatility
      %  quoted from Q1i(g)
        pd=makedist("normal","mu",0,"sigma",1);
        
%%
        d1=(log(V/D)+(rf+0.5*sigma^2)*T)/(sigma*sqrt(T));
        d2=d1-sigma*sqrt(T);
        Nd1=cdf(pd,d1);
        Nd2=cdf(pd,d2);
        NNd1=cdf(pd,-d1);
        NNd2=cdf(pd,-d2);
        
%%
        CallBS=V*Nd1-D*Nd2*exp(-rf*T);
        PutBS=D*NNd2*exp(-rf*T)-V*NNd1;
        % B_t=D*exp(-rf*T)-PutBS;
        B_t=V-CallBS;

%%          
            % SOLVE FOR CREDIT SPREAD [s]
            % D_t = K*exp(-(r+s)*(T-t))
            % s = -(log(D_t./K)./(T-t))-r;
            
        y=(log(D/B_t))/T;
        s=y-rf;
        results_a(i,j)=s;
        
    end
end
writematrix(results_a,'results_a.xlsx','Sheet',1);

%results_a =

%    0.3469    0.1779    0.0057    0.0000
%    0.1845    0.1096    0.0118    0.0002
%    0.1298    0.0833    0.0145    0.0007
%    0.1019    0.0687    0.0156    0.0013
%    0.0848    0.0592    0.0161    0.0020
%    0.0731    0.0524    0.0162    0.0026
%    0.0646    0.0473    0.0161    0.0031
%    0.0580    0.0433    0.0159    0.0036
%    0.0528    0.0400    0.0156    0.0039
%    0.0486    0.0373    0.0153    0.0042
%    0.0451    0.0350    0.0150    0.0045
%    0.0421    0.0330    0.0148    0.0047
%    0.0396    0.0312    0.0145    0.0049
%    0.0373    0.0297    0.0142    0.0051
%    0.0354    0.0284    0.0139    0.0052
%    0.0336    0.0271    0.0136    0.0053
%    0.0321    0.0260    0.0134    0.0054
%    0.0307    0.0250    0.0131    0.0055
%    0.0294    0.0241    0.0129    0.0055
%    0.0283    0.0233    0.0126    0.0056

% results=results*10000;

% plot V=80 and V=90
x=T_a;
y1=results_a(:,1);
y2=results_a(:,2);
xlabel('Maturity t(years)')
ylabel('Credit Spread (b.p.)')
plot(x,y1,x,y2)

% plot V=125 and V=180
x=T_a;
y3=results_a(:,3);
y4=results_a(:,4);
plot(x,y3,x,y4)

%%
%b

V_b=120;
sigma_b=[0.15,0.2,0.3,0.5];
T_b=0.5:0.5:10;
rf_b=0.028;
results_b=zeros(20,4); %the table output here

for i=1:20
    for j=1:4
      %  do the loop for matrix
        V=V_b;
        T=T_b(i);
        rf=rf_b;  %  rf=rfAll(j);
        D=90; %  face value of debt
        sigma=sigma_b(j); %volatility
      %  quoted from Q1i(g)
        pd=makedist("normal","mu",0,"sigma",1);
        d1=(log(V/D)+(rf+0.5*sigma^2)*T)/(sigma*sqrt(T));
        d2=d1-sigma*sqrt(T);
        Nd1=cdf(pd,d1);
        Nd2=cdf(pd,d2);
        NNd1=cdf(pd,-d1);
        NNd2=cdf(pd,-d2);
        CallBS=V*Nd1-D*Nd2*exp(-rf*T);
        PutBS=D*NNd2*exp(-rf*T)-V*NNd1;
        B_t=V-CallBS;            
        y=(log(D/B_t))/T;
        s=y-rf;
        results_b(i,j)=s;
        
    end
end
x=T_b;
y1=results_b(:,1);
y2=results_b(:,2);
y3=results_b(:,3);
y4=results_b(:,4);

plot(x,y1,x,y2,x,y3,x,y4)

%%
%c

results_c=zeros(20,2);
V_c=120;
sigma_c=0.18;
T_c=0.5:0.5:10;
rf_c=[0.03,0.09];

for i=1:20
    for j=1:2
        V=V_c;
        T=T_c(i);
        rf=rf_c(j);
        D=95;
        sigma=sigma_c;
        pd=makedist("normal","mu",0,"sigma",1);
        d1=(log(V/D)+(rf+0.5*sigma^2)*T)/(sigma*sqrt(T));
        d2=d1-sigma*sqrt(T);
        Nd1=cdf(pd,d1);
        Nd2=cdf(pd,d2);
        NNd1=cdf(pd,-d1);
        NNd2=cdf(pd,-d2);
        CallBS=V*Nd1-D*Nd2*exp(-rf*T);
        PutBS=D*exp(-rf*T)*NNd2-V*NNd1;
        B_t=V-CallBS;
        %B=D*exp(-rf*T)-PutBS;
        y=(log(D/B_t))/T;
        s=y-rf;
        results_c(i,j)=s;
    end
end
writematrix(results_c,'results_c.xlsx','Sheet',1);

% results=results*10000;
x=T_c;
y1=results_c(:,1);
y2=results_c(:,2);

plot(x,y1,x,y2)




%% SUBFUNCTIONS

function p=fcnN(x)
p=0.5*(1.+erf(x./sqrt(2)));
end
%
function p=fcnn(x)
p=exp(-0.5*x.^2)./sqrt(2*pi);
end

function Y = inv3d(X)
    Y = -X;
    Y(2,2,:) = X(1,1,:);
    Y(1,1,:) = X(2,2,:);
    detMat = 1./(X(1,1,:).*X(2,2,:) - X(1,2,:).*X(2,1,:));
    detMat = detMat(ones(1,2),ones(2,1),:);
    Y = detMat.*Y;
end
