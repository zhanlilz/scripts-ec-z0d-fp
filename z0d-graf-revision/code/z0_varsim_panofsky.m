function [z0,N,Nin] = z0_varsim_panofsky(data,z)%%%
%Estimate z from std(w) and u (plus u* and L for diagnosis)
% by flux-variance similarity, as in Panofsky 1984, Boundary-Layer Meteorol 28:305-308.
%   Input arguments:
%   data: m by 6 matrix of sigma(w)(=data(:,1)),u*(=data(:,2),L(=data(:,3),
%       u (=data(:,4))
%       measurements at an ECstation / sonic anemometer at m points in time
%   z: estimate of z-d


C1=1.3;%%%%C2=2.0;C3=0.99;C4=0.06; % estimates for the universal constants,
                % as used for this method by Toda & Sugita 2003
                % (Panofsky & Dutton 1984)
min_nobs = 10;

if size(data,2)~=4
    error('1st input argument must be m by 4 matrix')
end

flags=~isfinite(data);A=max(flags');data(A'==1,:)=[]; % only available data
% note that here u*, otherwise not used by the methods, plays an important
% role since the calling function (comparison of schemes to determine z0 / d)
% uses NaNning of u* values to perform the gerneral quality control (footprint etc.)

Nin=size(data,1);

flags=[...
      (~isfinite(data)) ... % listwise deletion of missing values
      (1./data(:,3)<-0.4)... % limits for L (must be near-neutral for Panofsky method)
      (1./data(:,3)>0.4)... 
      (data(:,4)<1.5)...
      ];
    
A=max(flags');

data(A'==1,:)=[];%,[],2),:)=[]; % compress matrix to used data

N=size(data,1);
if N<min_nobs %%%%
    warning('z0_varsim_panofsky: too few data points left after validation => NaN')
    z0=nan;
else

u=data(:,4);
sigw=data(:,1);
slope=u\sigw; % regression as in Eq. (2) of Panofsky, without offset as indicated by the figure
z0_z=1./exp(0.4*C1/slope);
z0=z0_z.*z;

end

