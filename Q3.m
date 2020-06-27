%% Author: Xinyu Zhang
%% CID: 01787342
%% Firm values in structured credit model
%
%

%%%%%%%%%%%%%%%%%%%%%%% LECTURE7 %%%%%%%%%%%%%%%%%%%%
SIGMA = [1.00 , 0.80 ; 0.80 , 1.00 ] ; 
% variance are 1.00 so off−diagonals are
% correlations % Our gaussian variables
R = mvnrnd(MU,SIGMA,1000);
% Run through S.N CDF
U = normcdf (R) ;
%Now U is uniform − check: 
figure ;
subplot (1 ,2 ,1)
histogram (U(: ,1) )
subplot (1 ,2 ,2) 
histogram (U(: ,2) )
% Now x is gamma distributed when we do this
x = gampdf(U(:,1),1,2);
% And y is beta distributed when we do this
y = betapdf(U(:,2),1,2);
% but they are correlated in the way we would like
figure ; 
scatter(x,y) 
%rho(x,y)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=250000;
corr=[0,0.15,0.35];

results_Q3=zeros(6,6);
for i=1:3
    results_Q3(:,i)=Q3_Gcopula(corr(i))';
end

for i=4:6
    results_Q3(:,i)=Q3_tcopula(corr(i-3))';
end


writematrix(results_Q3,'copula.xlsx','Sheet',1);







