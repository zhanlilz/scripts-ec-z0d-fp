function [z,z0, z2, z02, N, Nin] = z0d_mart(data, z0max, zmax, zv, min_nobs)%%%
%Estimate z=zm-d and z0 from timeseries of u, u*, and L - Martano type
%   Input arguments:
%   data: m by 3 matrix of u(=data(:,1)), u*(=data(:,2), L(=data(:,3)
%       measurements at an ECstation / sonic anemometer at m points in time
% version mart2, other than mart, gives back NaNs if the minimum was found
% at an edge of the screened range
% This version uses z0/L in addition to wind speed for quality control,
% for this an expected maximum z0 can be given as 2nd parameter

k=0.4;
% min_nobs = 10;

if size(data,2)~=3
    error('1st input argument must be m by 3 matrix')
end
if nargin<3
    zmax = 2.5; % 2.5 m
    if nargin <2
        z0max = 0.1; % 0.1 m
    end
end

flags=~isfinite(data);A=max(flags');data(A'==1,:)=[]; % only available data

Nin=size(data,1);

flags=[...
      (~isfinite(data)) ... % listwise deletion of missing values
      (data(:,1)<1.5) ... % lower wind speed boundary used by Martano
      (1./data(:,3)<-0.084./z0max) ...
      (1./data(:,3)>0.037./z0max) ...
      (1./data(:,3)>1./zmax) ...
      ];
    
A=max(flags');

data(A'==1,:)=[];%,[],2),:)=[]; % compress matrix to used data
N=size(data,1);
if N<min_nobs %%%%
    warning('z0d_mart: too few data points left after validation => NaN')
    z=nan;z0=nan;z2=nan;z02=nan;z3=nan;z03=nan;%%%
else

% mi=1;ma=60; % move zv to input argument with units in meters
% zv = (mi:ma)*0.1;
mi = 1; ma = numel(zv);

[zv_mat, lm_mat] = meshgrid(zv, data(:, 3));
[_, u_mat] = meshgrid(zv, data(:, 1));
[_, ustar_mat] = meshgrid(zv, data(:, 2)); 
[z0v, _, _, _, _, FI] = roughnesslength_revised3(zv_mat, lm_mat, u_mat, ustar_mat, k);
Sv = k .* u_mat ./ ustar_mat - FI;

z0avg=mean(z0v);
z0std=std(z0v);
z0cv=z0std./z0avg;
Svstd=std(Sv);
[dump,ix]=min(z0cv); 
z0=z0avg(ix);
z=zv(ix);
if ix==mi||ix==ma
    z=nan;
    z0=nan;
end

% Below corresponds to Eq.(7), FP-IT-1? even though it is named here z02 and
% z2? 
[dump,ix]=min(Svstd);
z02=z0avg(ix);
z2=zv(ix);
if ix==mi||ix==ma
    z2=nan;
    z02=nan;
end
    
end

