function [z0,phi,u,m,n,FI] = roughnesslength_revised3(zm,LM,U_meas,Ustar,k)

%%% 3rd revision: just for convenient use of the function within the
%%% Martano-type z0 d estimation => give out Fi also

%%% the coefficients of zm/L in this alternative rougness length calculator are changed
%%% from 5 (stable) and 16 (unstable) as given by Pauslon 1970
%%% to 6 and 19.3, respecively, in agreement with revisions to the Businger 1971 universal function
%%% by Högström 1988 (see e.g. Foken 2003, p. 44, or Table 2.8 in Foken 2008, p46)
%%% Foken, T., 2008. Micrometeorology. Springer Berlin Heidelberg, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-74666-9

  % another little revision, to account for z-less scaling
  sflag = zm ./ LM > 1;
  LM(sflag) = zm(sflag);

  z0 = zeros(size(zm));
  phi = zeros(size(zm));
  u = zeros(size(zm));
  m = zeros(size(zm));
  n = zeros(size(zm));
  FI = zeros(size(zm));

  sflag = LM > 0;
  FI(sflag) = 6 .* zm(sflag) ./ LM(sflag);
  phi(sflag) = 1 + 6 .* zm(sflag) ./ LM(sflag);
  alpha = U_meas(sflag) .* 0.4 ./ Ustar(sflag) - FI(sflag);
  z0(sflag) = zm(sflag) ./ exp(alpha);                           %function!!!
  %z0 = z0limit(z0);                               % set a lower limit for z0 to avoid numeric problems
  u(sflag) = Ustar(sflag) ./ k .* (log(zm(sflag) ./ z0(sflag)) + 5 .* zm(sflag) ./ LM(sflag));
  m(sflag) = (1 + 5 .* zm(sflag) ./ LM(sflag)) ./ (log(zm(sflag) ./ z0(sflag)) + 5 .* zm(sflag) ./ LM(sflag));
  n(sflag) = 1 ./ (1 + 5 .* zm(sflag) ./ LM(sflag));

  sflag = LM == 0;
  phi(sflag) = 1;
  alpha = U_meas(sflag) .* 0.4 ./ Ustar(sflag);
  z0(sflag) = zm(sflag) ./ exp(alpha);
  %z0 = z0limit(z0);
  u(sflag) = Ustar(sflag) ./ k .* (log(zm(sflag) ./ z0(sflag)));
  m(sflag) = Ustar(sflag) ./ k ./ u(sflag);
  n(sflag) = 1;

  sflag = LM < 0;
  zeta = (1 - (19.3 .* zm(sflag) ./ LM(sflag))).^0.25;
  FI(sflag) = -2 .* log((1 + zeta) / 2) - log((1 + zeta .* zeta) / 2) + 2 .* atan(zeta) - pi / 2;
  phi(sflag) = 1 ./ zeta.^ 2;
  alpha = U_meas(sflag) * 0.4 ./ Ustar(sflag) - FI(sflag);
  z0(sflag) = zm(sflag) ./ exp(alpha);
  %z0 = z0limit(z0);
  u(sflag) = Ustar(sflag) ./ k .* (log(zm(sflag) ./ z0(sflag)) + FI(sflag));
  m(sflag) = Ustar(sflag) ./ k ./ zeta ./ u(sflag);
  n(sflag) = (1 - 24 .* zm(sflag) ./ LM(sflag)) ./ (1 - 16 .* zm(sflag) ./ LM(sflag));
