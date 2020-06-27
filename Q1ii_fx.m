
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

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


