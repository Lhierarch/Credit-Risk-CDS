%% Author: Xinyu Zhang
%% CID: 01787342
%
%

global C;
    global t;
    global S;
%GM C=[1.60181,1.29915,1.13627,1.02703,0.9575,0.910195,0.86289,0.837007,0.8111233,0.78524];
%Ford 
C=[0.90505,0.75648,0.69219,0.64886,0.60888,0.587895,0.56691,0.5577967,0.5486833,0.53957];
%DB C=[0.01174,0.01277,0.0132,0.01308,0.0130,0.012985,0.01297,0.0129533,0.0129367,0.01292];
%Kraft C=[0.00633,0.00726,0.00825,0.00914,0.00999,0.01036,0.01073,0.010913,0.011097,0.01128];
%GMAC C=[0.22516,0.16252,0.13703,0.12298,0.11364,0.11072,0.1078,0.10354,0.09928,0.09502];
%JPM C=[0.01075,0.01134,0.01161,0.01204,0.0125,0.01275,0.0130,0.0130133,0.0130267,0.01304];
%Heinz C=[0.00526,0.00628,0.00732,0.00838,0.00939,0.009605,0.00982,0.0102,0.01058,0.01096];


   
%%
%a
%pick GM
%C=[1.60181,1.29915,1.13627,1.02703,0.9575,0.910195,0.86289,0.837007,0.8111233,0.78524];

options = optimset('PlotFcns','optimplotfval','TolFun',1e-7);
fun=@Q2_fx; % my own function
x0=[1];
[x,fval] = fminsearch(fun,x0,options);
S=zeros(1,11);
S(1)=1;
f=zeros(1,11);

for t=1:10
    [x,fval]=fminsearch(fun,x0,options);
    S(t+1)=x;
    f(t+1)=fval;
end


lambda=zeros(1,11);
conditionalPD=zeros(1,11);

for i =1:10
   lambda(i+1)=-log(S(i+1)/S(i));
   % Conditional PD=exp⁡(-λ_t*1)
   conditionalPD(i+1)=exp(-lambda(i+1)*1);
end

ST=transpose(S);
LT=transpose(lambda);
CT=transpose(conditionalPD)

%writematrix(ST,'ST.xlsx','Sheet',1);
%writematrix(LT,'ST.xlsx','Sheet',2);
writematrix(CT,'CT.xlsx','Sheet',1);


