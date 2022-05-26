function [z1,z2,N,Nin] = z0d_varsim(data, zmax, zv, min_nobs)%%%
%Estimate z or d from std(w), u*, and L - variance similarity type
%   Input arguments:
%   data: m by 6 matrix of sigma(w)(=data(:,1)),u*(=data(:,2),L(=data(:,3),
%       u (=data(:,4)), sigma(T) (=data(:,5)), T* (=data(:,6))
%       measurements at an ECstation / sonic anemometer at m points in time

C1=1.3;C2=2.0;C3=0.99;C4=0.06; % estimates for the universal constants,
                % as used for this method by Toda & Sugita 2003
                % (Panofsky & Dutton 1984)
% min_nobs = 10;

if size(data,2)~=6
    error('1st input argument must be m by 6 matrix')
end
if nargin<2
    zmax = 2.5; % 2.5 m
end

flags=~isfinite(data);A=max(flags');data(A'==1,:)=[]; % only available data

Nin=size(data,1);

flags=[...
      (~isfinite(data)) ... % listwise deletion of missing values
      (zmax./data(:,3)>0)... % limits for L
      (zmax./data(:,3)<-99999)...
      (data(:,2)<0.05)... % further limits for u* and t* found necessary
      (abs(data(:,6))<0.3)...% by ourselves to obtain robust results
      ];  
    
A=max(flags');

data(A'==1,:)=[];%,[],2),:)=[]; % compress matrix to used data

N=size(data,1);
if N<min_nobs %%%%
    warning('z0d_varsim: too few data points left after validation => NaN')
    z1=nan;z2=nan;%%%
else

% mi=1;ma=60;
% zv = (mi:ma)*0.1;
mi = 1; ma = numel(zv);

[zv_mat, lm_mat] = meshgrid(zv, data(:, 3));
[_, stdw_mat] = meshgrid(zv, data(:, 1));
[_, ustar_mat] = meshgrid(zv, data(:, 2)); 
[_, stdt_mat] = meshgrid(zv, data(:, 5));
[_, tstar_mat] = meshgrid(zv, data(:, 6));

lhs1 = stdw_mat ./ ustar_mat; %Left hand side for vertical wind relation
lhs2 = stdt_mat ./ abs(tstar_mat); % ...for (virtual) temperature relation

rhs1 = C1 .* nthroot((1 - C2 .* zv_mat ./ -abs(lm_mat)), 3); % right hand sides...
rhs2 = C3 .* (C4-zv_mat ./ -abs(lm_mat)).^(-1 / 3);
SE1 = (rhs1 - lhs1).^2;
SE2 = ((rhs2 - lhs2) ./ rhs2).^2;

MSE=(mean(SE1));
%figure;
%plot(zv,MSE,'+');
[dump,ix]=min(MSE);
z1=zv(ix);
if ix==mi||ix==ma
    z1=nan;
end

MSE=(mean(SE2));
%figure;
%plot(zv,MSE,'+');
[dump,ix]=min(MSE);
z2=zv(ix);
if ix==mi||ix==ma
    z2=nan;
end
   
end

