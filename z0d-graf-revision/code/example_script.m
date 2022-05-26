% Matlab/Octave script for simultaneous estimation of roughness length z0
% and effective measurement height z (and thus d, if zm is known)
% from sonic anemometer eddy covariance measurements of friction velocity,
% mean wind speed, and buoyancy flux at a single level but over a range
% of stabilities. The actual algorithms are in functions starting z0d_...
% A detailed description of what´s happening is given in:
% Graf et al., Intercomparison of methods for the simultaneous estimation
% of displacement height and roughness length from single-level eddy
% covariance data. Boundary-Layer Meteorol., in press.
% a.graf@fz-juelich.de, January 2011 - September 2013.

clear; clc(); close all;

% User-changeable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input file(s)
path='../test-inputs/ECdata_'; % path and prefix
datafiles={...
            'J1lower'
            }; %names of datafile(s) without file extension
ext='.txt'; % file extension
sep=',';  % column separator
header=1;   % number of header lines to omit
miss=-9999; % entry that indicates "missing data" (apart from NaN)
% which column contains...?
col_ind=1; % any independent variable that might affect z and z0*
col_u=7;    % the mean wind speed (m/s)
col_ustar=9;% friction velocity (m/s)
col_Hv=6;   % the buoyancy flux (W/m²) (0 if not separately available)
col_Tv=4;   % the virtual temperature (K) (0 if not separately available)
col_rho=5;  % the mean air density (kg/m³) (0 if not separately available)
col_L=0;    % Obukhov length L (Default 0, only give if Hv,Tv,rho missing)
col_stdw=3; % needed for the flux-variance similarity method
col_stdT=10; % ".
% *ind (with indstep) will control the definition of sub-datasets within
% each file. Results will be plotted as a function of ind. ind can be:
% - the Day of Year (e.g. col 1), if the canopy changes throughout the year
% - the geometric measurement height above ground, if it changes
% - Canpoy height or density, if such measurements are available
% - the expected effective height (e.g. col13), as a result of both 
% - wind direction (e.g. col12), if the canopy differs in each direction
% To calculate only one set of results for the whole dataset, set col_ind=0
%
% "constants":
indstep=30; % binning width for ind. Enlarge if results scatter too much
running=1;  % 1 for running window across ind of width indstep, 0 for block
breakpoint=[216];
    % 1 point in 'ind'allowed for each station (input file),
    % where data below and above this point will never be in the same bin.
    % e.g., enter day of harvest, if col_ind points to the DOY column.
    % Enter anything above max(ind) to deactivate the breakpoint.
mergepoint=[999]; %similarly, define one point
    % at the upper end of ind, that will be interpreted as 0 if necessary.
    % E.g. if ind is wind direction, 360°. 
    % Enter anything > max(ind)+indstep to deactivate (e.g. 999).
    % mergepoint won´t work if breakpoint is used!
col_exclude=8; % column used to exclude invalid data (e.g. wind direction
            % column to exclude the sector of the mast / anemometer back).
lower=[0]; % lower limit (for each station) of col_exclude
upper=[360];               % upper
col_exclude2=11; % 2nd exclusion, e.g. footprint target contribution
lower2=[0.5]; %
upper2=[999]; %
col_exclude3=12; % in this example, combined footprint contribution
lower3=[0.9]; % of low crop heights (sugarbeet + wheat + barley)
upper3=[999];
Tv=NaN;     % only needed if neither Tv nor L given in input file
rho=NaN;    % same. Assuming constant values may cause systematic L errors!
g=9.81;     % gravity acceleration (9.81 m/s²)
k=0.4;      % von Karman constant (0.4)
cp=1005;    % specific heat capacity of air @ const. pressure (1005 J/kg/K)

% Readin & homogenize part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hlp=floor(log10(indstep)); % help variable to group similar 'ind' below the
hlp=10^(-hlp+1); % ...indstep-level for speed and figure clarity.
nos=size(datafiles,1); %number of stations/files
for j=1:nos % loop through stations
    data=dlmread(strcat(path,char(datafiles(j)),ext),sep,header,0);
    recs=size(data,1); % number of records
    data(data==miss)=NaN;
    exclude=data(:,col_exclude)<lower(j)|data(:,col_exclude)>upper(j);
    data(exclude,col_ustar)=NaN;
    exclude2=data(:,col_exclude2)<lower2(j)|data(:,col_exclude2)>upper2(j);
    data(exclude2,col_ustar)=NaN;
    exclude3=data(:,col_exclude3)<lower3(j)|data(:,col_exclude3)>upper3(j);
    data(exclude3,col_ustar)=NaN;
    clear u ustar L zL indbin x y reslts;
    u=data(:,col_u);
    ustar=data(:,col_ustar);
    if col_Hv>0; clear Hv; Hv=data(:,col_Hv); end
    if col_Tv>0; clear Tv; Tv=data(:,col_Tv);
        else Tv=repmat(Tv,recs,1); end
    if col_rho>0; clear rho; rho=data(:,col_rho);
        else rho=repmat(rho,recs,1); end
    if col_ind>0; ind=data(:,col_ind); else ind=ones(recs,1); end
    ind=round(ind*hlp)/hlp; % see remarks on 'round'
    if col_L>0; L=data(:,col_L); else L=NaN(recs,1); end
    Tstar=NaN(recs,1); 
    stdw=data(:,col_stdw);
    stdT=data(:,col_stdT);
    for i=1:recs % loop through timestamps
        if isfinite(Hv(i)) && isfinite(Tv(i)) ...
                && isfinite(rho(i)) && isfinite(ustar(i))
            % compute Obukhov length
            L(i)=-rho(i)*(ustar(i))^3/(k*g*((Hv(i))/(cp*Tv(i))));
            % compute T*
            Tstar(i) = - Hv(i)/ ( rho(i) * cp * ustar(i) );
        end
        if running==0
            indbin(i)=int16(ind(i)/indstep); %prepare division into bins
        else
            indbin(i)=int16(ind(i)*hlp); % or prepare for running window
        end
    end
    if running==0
        tol=0;
        mergepoint(j)=mergepoint(j)/indstep;
    else
        tol=indstep*hlp;
        mergepoint(j)=mergepoint(j)*hlp;
    end
    indbin=indbin-min(indbin)+1;
    mergepoint(j)=mergepoint(j)-min(indbin)+1;
    indbin(ind>breakpoint(j))=indbin(ind>breakpoint(j))+tol/2+1;
 
% Call analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for h=1:max(indbin) % loop through the bins
        flags=((abs(indbin-h)<=tol/2|...
            abs(indbin-mergepoint(j)-h)<=tol/2|...
            abs(indbin+mergepoint(j)-h)<=tol/2));
        u_=u(flags);
        ustar_=ustar(flags);
        L_=L(flags);
        stdw_=stdw(flags);
        stdT_=stdT(flags);
        Tstar_=Tstar(flags);
        Ntotal(h)=size(u_,1);
        % 1st model set - linear simple and 3d regression  
        [z1(h),z01(h),z2(h),z02(h),Nreg(h),Ninreg(h)]= z0d_reg([u_ ustar_ L_], 0.1, 2.5, 30);
        % 2nd set of models - iterative Martano   
        [z3(h),z03(h), z4(h), z04(h), Nmart(h),Ninmart(h)]=z0d_mart([u_ ustar_ L_], 0.1, 2.5, (0.1:0.1:6), 30);
        % 3rd set of models - Weaver/Rotach/TodaSugita type (only d)
        [z5(h),z6(h),Nvar(h),Ninvar(h)]=z0d_varsim([stdw_ ustar_ L_ u_ stdT_ Tstar_], 2.5, (0.1:0.1:6), 30);
        % estimate z0 from variance(w) (after estimating d above)
        [z05(h),Nvarz0(h),Ninvarz0(h)]=z0_varsim_panofsky([stdw_ ustar_ L_ u_], z5(h));
        if ( Ninreg(h)~=Ninmart(h) || Ninreg(h)~=Ninvar(h) )
            [Ninreg(h) Ninmart(h) Ninvar(h)]
        end
    end;

% Output part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:recs
        reslts(i,:)=[z1(indbin(i)) z01(indbin(i)) ...
            z2(indbin(i)) z02(indbin(i))...
            z3(indbin(i)) z03(indbin(i)) ...
            z4(indbin(i)) z04(indbin(i)) ...
            z5(indbin(i)) z05(indbin(i)) ...
            z6(indbin(i)) ...
            Nreg(indbin(i)) Nmart(indbin(i)) Nvar(indbin(i)) Nvarz0(indbin(i)) Ntotal(indbin(i)) Ninreg(indbin(i))];
    end
    data=[data reslts];
    outfname=strcat('../test-outputs/',char(datafiles(j)),'_withz0d','.csv');
    fid=fopen(outfname, 'wt');
    fprintf(fid, ['DOY,hourmin,std(w),Tv,Rhoson,Hv,u,dir,u*,std(T),FP_target,FP_all,' ...
                  'z1,z01,z2,z02,z3,z03,z4,z04,z5,z05,z6,Nreg,Nmart,Nvar,Nvarz0,Ntotal,Ninreg\n'])
    fclose(fid);
    dlmwrite(outfname,data,'-append','precision','%.6f');
    subplot(3,nos,j);
        axis off;
        text(.5,.5,{char(datafiles(j))},'HorizontalAlignment', 'center');
    subplot(3,nos,nos+j);
        plot(ind,reslts(:,1),'.'...
            ,ind,reslts(:,3),'+'...
            ,ind,reslts(:,5),'*'...
            ,ind,reslts(:,7),'o'...
            ,ind,reslts(:,9),'-'...
            ,ind,reslts(:,10),'-'...
            );
    subplot(3,nos,2*nos+j);
        semilogy(ind,reslts(:,2),'.'...
            ,ind,reslts(:,4),'+'...
            ,ind,reslts(:,6),'*'...
            ,ind,reslts(:,8),'o');
end
    


