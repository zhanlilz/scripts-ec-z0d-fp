function [z1,z01,z2,z02,N,Nin] = z0d_reg(data,z0max,zmax,min_nobs)
%Estimate z=zm-d and z0 from timeseries of u, u*, and L - regression
%   Input arguments:
%   data: m by 3 matrix of u(=data(:,1)), u*(=data(:,2), L(=data(:,3)
%       measurements at an ECstation / sonic anemometer at m points in time
%   bounds (optional): 2 element vector of lower & upper bound for 1/L.
%
%   Despcription
%   Starting from the logarithmic wind profile
%   u=u*/k*(ln(z/z0)-f(z/L)) where
%       u wind speed m/s
%       u* friction velocity m/s
%       k von karman constant
%       z measurement height (minus any displacement height) m
%       z0 roughness length
%       L Obukhoc length m
%       f Intergal of the universal function for momentum
%   and confining to the near-neutral to moderately stable range, i.e.
%       L > max(bounds) or L < min(bounds) with
%       max(bounds) > zest/(0.5...1)
%       min(bounds) < zest/(-0.05...0)
%       zest = conservative estimate (upper limit) to z,
%   then f(z/L)~-uni*(z/L) (uni = 6 according to current literature)
%   and we can put up a regression
%       uk/u*=a+b(1/L)
%   to determine z (and from that d) as b/uni and z0 as z/exp(a)
%   a.graf@fz-juelich.de

k=0.4;uni=6;
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
    (1./data(:,3)<-0.103./zmax) ... 
    (1./data(:,3)<-0.084./z0max) ... % check stability range
    (1./data(:,3)>0.037./z0max) ...
    (1./data(:,3)>1./zmax) ...
    (~isfinite(data)) ... % listwise deletion of missing values
    ];
    
A=max(flags');

data(A'==1,:)=[];%,[],2),:)=[]; % compress matrix to used data
N=size(data,1);
if N<min_nobs %%%%
    warning('z0d_reg: too few data points left after validation => NaN')
    z1=nan;z01=nan;z2=nan;z02=nan;%%z3=nan;z03=nan;z4=nan;z04=nan;
else
        
    % 1st model (linear regression)
    y=data(:,1).*k./data(:,2);  % u*k/u*
    x=uni./data(:,3);       % uni/L
    B=[ones(size(x)) x]\y;
    z1=B(2);
    z01=z1/exp(B(1));
    
    % 2nd model (3D regression without offset)
    clear B; %%y1 x1 x2;
    y1=data(:,1); % u
    x1=data(:,2)./k; % u*/k
    x2=x.*x1; % uni*u*/(k*L), or uni*g*Hv/(Tv*rhocp*u*²)
    B=[x1 x2]\y1;
    z2=B(2);
    z02=z2/exp(B(1));
          
end

