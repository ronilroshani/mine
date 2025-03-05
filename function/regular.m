function [out]=regular(parameter1,parameter2,method,action)
% parameter1=design matrix
% parameter2=observation
% method: 'GCV'-'LC'
% action: 'tikh','tsvd','dsvd'
%% ------------------------------------------------------------------------

s=svd(parameter1);
[U,~,VV]=svd(parameter1);
if strcmp(method,'GCV')
[reg_min,~,~] = gcv(U,s,parameter2,action);
[out,~,~] = tikhonov(U,s,VV,parameter2,reg_min);
elseif strcmp(method,'LC')
[reg_corner,~,~,~]=l_curve(U,s,parameter2,action);
[out,~,~] = tikhonov(U,s,VV,parameter2,reg_corner);
end

function [reg_min,G,reg_param] = gcv(U,s,b,method)
%GCV Plot the GCV function and find its minimum.
%
% [reg_min,G,reg_param] = gcv(U,s,b,method)
% [reg_min,G,reg_param] = gcv(U,sm,b,method)  ,  sm = [sigma,mu]
%
% Plots the GCV-function
%          || A*x - b ||^2
%    G = -------------------
%        (trace(I - A*A_I)^2
% as a function of the regularization parameter reg_param.
% Here, A_I is a matrix which produces the regularized solution.
%
% The following methods are allowed:
%    method = 'Tikh' : Tikhonov regularization   (solid line )
%    method = 'tsvd' : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd' : damped SVD or GSVD        (dotted line)
% If method is not specified, 'Tikh' is default.
%
% If any output arguments are specified, then the minimum of G is
% identified and the corresponding reg. parameter reg_min is returned.


%% ------------------------------------------------------------------------

% Set defaults.
if (nargin==3), method='Tikh'; end  % Default method.
npoints = 200;                      % Number of points on the curve.
smin_ratio = 16*eps;                % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(s);
beta = U'*b; beta2 = norm(b)^2 - norm(beta)^2;
if (ps==2)
  s = s(p:-1:1,1)./s(p:-1:1,2); beta = beta(p:-1:1);
end
if (nargout > 0), find_min = 1; else find_min = 0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!!!!!!!!!!!!!!!!!!
    GGG=@gcvfun;   %baraye tarif gcv in khat ezafe shod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4)) 

   
  % Vector of regularization parameters.
  reg_param = zeros(npoints,1); G = reg_param; s2 = s.^2;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Intrinsic residual.
  delta0 = 0;
  if (m > n & beta2 > 0), delta0 = beta2; end

  % Vector of GCV-function values.
  for i=1:npoints
    G(i) = gcvfun(reg_param(i),s2,beta(1:p),delta0,m-n);
  end 

  % Plot GCV function.
%   figure(1)
%   loglog(reg_param,G,'-'), xlabel('\lambda'), ylabel('G(\lambda)')
%   title('GCV function')

  % Find minimum, if requested.
  if (find_min)
    [minG,minGi] = min(G); % Initial guess.
    reg_min = fminbnd(GGG,...
      reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
      optimset('Display','off'),s2,beta(1:p),delta0,m-n); % Minimizer.
    minG = gcvfun(reg_min,s2,beta(1:p),delta0,m-n); % Minimum of GCV function.
%     ax = axis;
%     HoldState = ishold; hold on;
%     loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
% %     title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
%     axis(ax)
%     if (~HoldState), hold off; end
  end

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
   
  % Vector of GCV-function values.
  rho2(p-1) = abs(beta(p))^2;
  if (m > n & beta2 > 0), rho2(p-1) = rho2(p-1) + beta2; end
  for k=p-2:-1:1, rho2(k) = rho2(k+1) + abs(beta(k+1))^2; end
  G = zeros(p-1,1);
  for k=1:p-1
    G(k) = rho2(k)/(m - k + (n - p))^2;
  end
  reg_param = (1:p-1)';

  % Plot GCV function.
%   semilogy(reg_param,G,'o'), xlabel('k'), ylabel('G(k)')
%   title('GCV function')

  % Find minimum, if requested.
  if (find_min)
    [minG,reg_min] = min(G);
%     ax = axis;
%     HoldState = ishold; hold on;
%     semilogy(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
% %     title(['GCV function, minimum at k = ',num2str(reg_min)])
%     axis(ax);
%     if (~HoldState), hold off; end
  end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  % Vector of regularization parameters.
  reg_param = zeros(npoints,1); G = reg_param;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Intrinsic residual.
  delta0 = 0;
  if (m > n & beta2 > 0), delta0 = beta2; end

  % Vector of GCV-function values.
  for i=1:npoints
    G(i) = gcvfun(reg_param(i),s,beta(1:p),delta0,m-n,1);
  end

  % Plot GCV function.
%   loglog(reg_param,G,':'), xlabel('\lambda'), ylabel('G(\lambda)')
%   title('GCV function')

  % Find minimum, if requested.
  if (find_min)
    [minG,minGi] = min(G); % Initial guess.
    reg_min = fminbnd(GGG,...
      reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
      optimset('Display','off'),s,beta(1:p),delta0,m-n,1); % Minimizer.
    minG = gcvfun(reg_min,s,beta(1:p),delta0,m-n,1); % Minimum of GCV function.
%     ax = axis;
%     HoldState = ishold; hold on;
%     loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
% %     title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
%     axis(ax)
%     if (~HoldState), hold off; end
  end

  

  
  
elseif (strncmp(method,'mtsv',4) | strncmp(method,'ttls',4))

  error('The MTSVD and TTLS methods are not supported')

else
  error('Illegal method')
  close all
end
function G = gcvfun(lambda,s2,beta,delta0,mn,dsvd)

% Auxiliary routine for gcv.  PCH, IMM, Feb. 24, 2008.

% Note: f = 1 - filter-factors.
if (nargin==5)
   f = (lambda^2)./(s2 + lambda^2);
else
   f = lambda./(s2 + lambda);
end
G = (norm(f.*beta)^2 + delta0)/(mn + sum(f))^2;
end
end
function [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
%TIKHONOV Tikhonov regularization.
%
% [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
% [x_lambda,rho,eta] = tikhonov(U,sm,X,b,lambda,x_0) ,  sm = [sigma,mu]
%
% Computes the Tikhonov regularized solution x_lambda.  If the SVD
% is used, i.e. if U, s, and V are specified, then standard-form
% regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || x - x_0 ||^2 } .
% If, on the other hand, the GSVD is used, i.e. if U, sm, and X are
% specified, then general-form regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || L (x - x_0) ||^2 } .
%
% If x_0 is not specified, then x_0 = 0 is used
%
% Note that x_0 cannot be used if A is underdetermined and L ~= I.
%
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Per Christian Hansen, IMM, April 14, 2003.

% Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of
% Ill-Posed Problems", Wiley, 1977.

% Initialization.
if (min(lambda)<0)
  error('Illegal regularization parameter lambda')
end
m = size(U,1);
n = size(V,1);
[p,ps] = size(s);
beta = U(:,1:p)'*b;
zeta = s(:,1).*beta;
ll = length(lambda); x_lambda = zeros(n,ll);
rho = zeros(ll,1); eta = zeros(ll,1);

% Treat each lambda separately.
if (ps==1)

  % The standard-form case.
  if (nargin==6), omega = V'*x_0; end
  for i=1:ll
    if (nargin==5)
      x_lambda(:,i) = V(:,1:p)*(zeta./(s.^2 + lambda(i)^2));
      rho(i) = lambda(i)^2*norm(beta./(s.^2 + lambda(i)^2));
    else
      x_lambda(:,i) = V(:,1:p)*...
        ((zeta + lambda(i)^2*omega)./(s.^2 + lambda(i)^2));
      rho(i) = lambda(i)^2*norm((beta - s.*omega)./(s.^2 + lambda(i)^2));
    end
    eta(i) = norm(x_lambda(:,i));
  end
  if (nargout > 1 & size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

elseif (m>=n)

  % The overdetermined or square general-form case.
  gamma2 = (s(:,1)./s(:,2)).^2;
  if (nargin==6), omega = V\x_0; omega = omega(1:p); end
  if (p==n)
    x0 = zeros(n,1);
  else
    x0 = V(:,p+1:n)*U(:,p+1:n)'*b;
  end
  for i=1:ll
    if (nargin==5)
      xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
      x_lambda(:,i) = V(:,1:p)*xi + x0;
      rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
    else
      xi = (zeta + lambda(i)^2*(s(:,2).^2).*omega)./...
           (s(:,1).^2 + lambda(i)^2*s(:,2).^2);
      x_lambda(:,i) = V(:,1:p)*xi + x0;
      rho(i) = lambda(i)^2*norm((beta - s(:,1).*omega)./...
               (gamma2 + lambda(i)^2));
    end
    eta(i) = norm(s(:,2).*xi);
  end
  if (nargout > 1 & size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

else

  % The underdetermined general-form case.
  gamma2 = (s(:,1)./s(:,2)).^2;
  if (nargin==6), error('x_0 not allowed'), end
  if (p==m)
    x0 = zeros(n,1);
  else
    x0 = V(:,p+1:m)*U(:,p+1:m)'*b;
  end
  for i=1:ll
    xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
    eta(i) = norm(s(:,2).*xi);
  end

end
end
function [reg_corner,rho,eta,reg_param] = l_curve(U,sm,b,method,L,V)
%L_CURVE Plot the L-curve and find its "corner".
%
% [reg_corner,rho,eta,reg_param] =
%                  l_curve(U,s,b,method)
%                  l_curve(U,sm,b,method)  ,  sm = [sigma,mu]
%                  l_curve(U,s,b,method,L,V)
%
% Plots the L-shaped curve of eta, the solution norm || x || or
% semi-norm || L x ||, as a function of rho, the residual norm
% || A x - b ||, for the following methods:
%    method = 'Tikh'  : Tikhonov regularization   (solid line )
%    method = 'tsvd'  : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd'  : damped SVD or GSVD        (dotted line)
%    method = 'mtsvd' : modified TSVD             (x markers  )
% The corresponding reg. parameters are returned in reg_param.  If no
% method is specified then 'Tikh' is default.  For other methods use plot_lc.
%
% Note that 'Tikh', 'tsvd' and 'dsvd' require either U and s (standard-
% form regularization) or U and sm (general-form regularization), while
% 'mtvsd' requires U and s as well as L and V.
%
% If any output arguments are specified, then the corner of the L-curve
% is identified and the corresponding reg. parameter reg_corner is
% returned.  Use routine l_corner if an upper bound on eta is required.

% Reference: P. C. Hansen & D. P. O'Leary, "The use of the L-curve in
% the regularization of discrete ill-posed problems",  SIAM J. Sci.
% Comput. 14 (1993), pp. 1487-1503.

% Per Christian Hansen, IMM, July 26, 2007.

% Set defaults.
if (nargin==3), method='Tikh'; end  % Tikhonov reg. is default.
npoints = 200;  % Number of points on the L-curve for Tikh and dsvd.
smin_ratio = 16*eps;  % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(sm);
if (nargout > 0), locate = 1; else locate = 0; end
beta = U'*b; beta2 = norm(b)^2 - norm(beta)^2;
if (ps==1)
  s = sm; beta = beta(1:p);
else
  s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
end

xi = beta(1:p)./s;


if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))

  eta = zeros(npoints,1); rho = eta; reg_param = eta; s2 = s.^2;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  for i=1:npoints
    f = s2./(s2 + reg_param(i)^2);
    eta(i) = norm(f.*xi);
    rho(i) = norm((1-f).*beta(1:p));
  end
  if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
  marker = '-'; txt = 'Tikh.';

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))

  eta = zeros(p,1); rho = eta;
  eta(1) = abs(xi(1))^2;
  for k=2:p, eta(k) = eta(k-1) + abs(xi(k))^2; end
  eta = sqrt(eta);
  if (m > n)
    if (beta2 > 0), rho(p) = beta2; else rho(p) = eps^2; end
  else
    rho(p) = eps^2;
  end
  for k=p-1:-1:1, rho(k) = rho(k+1) + abs(beta(k+1))^2; end
  rho = sqrt(rho);
  reg_param = (1:p)'; marker = 'o';
  if (ps==1)
    U = U(:,1:p); txt = 'TSVD';
  else
    U = U(:,1:p); txt = 'TGSVD';
  end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  eta = zeros(npoints,1); rho = eta; reg_param = eta;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  for i=1:npoints
    f = s./(s + reg_param(i));
    eta(i) = norm(f.*xi);
    rho(i) = norm((1-f).*beta(1:p));
  end
  if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
  marker = ':';
  if (ps==1), txt = 'DSVD'; else txt = 'DGSVD'; end

elseif (strncmp(method,'mtsv',4))

  if (nargin~=6)
    error('The matrices L and V must also be specified')
  end
  [p,n] = size(L); rho = zeros(p,1); eta = rho;
  [Q,R] = qr(L*V(:,n:-1:n-p),0);
  for i=1:p
    k = n-p+i;
    Lxk = L*V(:,1:k)*xi(1:k);
    zk = R(1:n-k,1:n-k)\(Q(:,1:n-k)'*Lxk); zk = zk(n-k:-1:1);
    eta(i) = norm(Q(:,n-k+1:p)'*Lxk);
    if (i < p)
      rho(i) = norm(beta(k+1:n) + s(k+1:n).*zk);
    else
      rho(i) = eps;
    end
  end
  if (m > n & beta2 > 0), rho = sqrt(rho.^2 + beta2); end
  reg_param = (n-p+1:n)'; txt = 'MTSVD';
  U = U(:,reg_param); sm = sm(reg_param);
  marker = 'x'; ps = 2;  % General form regularization.

else
  error('Illegal method')
end

% Locate the "corner" of the L-curve, if required.
if (locate)
  [reg_corner,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,sm,b,method);
end

% Make plot.
plot_lc(rho,eta,marker,ps,reg_param);
if locate
  ax = axis;
  HoldState = ishold; hold on;
  loglog([min(rho)/100,rho_c],[eta_c,eta_c],':r',...
         [rho_c,rho_c],[min(eta)/100,eta_c],':r')
  title(['L-curve, ',txt,' corner at ',num2str(reg_corner)]);
  axis(ax)
  if (~HoldState), hold off; end
end


function [reg_c,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,s,b,method,M)
%L_CORNER Locate the "corner" of the L-curve.
%
% [reg_c,rho_c,eta_c] =
%        l_corner(rho,eta,reg_param)
%        l_corner(rho,eta,reg_param,U,s,b,method,M)
%        l_corner(rho,eta,reg_param,U,sm,b,method,M) ,  sm = [sigma,mu]
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively (such as the output from routine l_curve).
%
% If nargin = 3, then no particular method is assumed, and if
% nargin = 2 then it is issumed that reg_param = 1:length(rho).
%
% If nargin >= 6, then the following methods are allowed:
%    method = 'Tikh'  : Tikhonov regularization
%    method = 'tsvd'  : truncated SVD or GSVD
%    method = 'dsvd'  : damped SVD or GSVD
%    method = 'mtsvd' : modified TSVD,
% and if no method is specified, 'Tikh' is default.  If the Spline Toolbox
% is not available, then only 'Tikh' and 'dsvd' can be used.
%
% An eighth argument M specifies an upper bound for eta, below which
% the corner should be found.

% Per Christian Hansen, IMM, July 26, 2007.

% Set default regularization method.
if (nargin <= 3)
  method = 'none';
  if (nargin==2), reg_param = (1:length(rho))'; end
else
  if (nargin==6), method = 'Tikh'; end
end

% Set this logical variable to 1 (true) if the corner algorithm
% should always be used, even if the Spline Toolbox is available.
alwayscorner = 0;

% Set threshold for skipping very small singular values in the
% analysis of a discrete L-curve.
s_thr = eps;  % Neglect singular values less than s_thr.

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.

% Initialization.
if (length(rho) < order)
  error('Too few data points for L-curve analysis')
end
if (nargin > 3)
  [p,ps] = size(s); [m,n] = size(U);
  size(U)
  size(b)
  beta = U'*b;
  size(beta)
  if (m>n), b0 = b - U*beta; end
  if (ps==2)
    s = s(p:-1:1,1)./s(p:-1:1,2);
    beta = beta(p:-1:1);
  end
  beta = beta(1:p);
  xi = beta./s;
end

% Restrict the analysis of the L-curve according to M (if specified).
if (nargin==8)
  index = find(eta < M);
  rho = rho(index); eta = eta(index); reg_param = reg_param(index);
end

if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))

  % The L-curve is differentiable; computation of curvature in
  % log-log scale is easy.

  % Compute g = - curvature of L-curve.
  g = lcfun(reg_param,s,beta,xi);
GGG=@lcfun;
  % Locate the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  [gmin,gi] = min(g);
  reg_c = fminbnd(GGG,...
    reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
    optimset('Display','off'),s,beta,xi); % Minimizer.
  kappa_max = - lcfun(reg_c,s,beta,xi); % Maximum curvature.

  if (kappa_max < 0)
    lr = length(rho);
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    f = (s.^2)./(s.^2 + reg_c^2);
    eta_c = norm(f.*xi);
    rho_c = norm((1-f).*beta);
    if (m>n), rho_c = sqrt(rho_c^2 + norm(b0)^2); end
  end

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4) | ...
        strncmp(method,'mtsv',4) | strncmp(method,'none',4))

  % Use the adaptive pruning algorithm to find the corner, if the
  % Spline Toolbox is not available.
  if ~exist('splines','dir') | alwayscorner
    %error('The Spline Toolbox in not available so l_corner cannot be used')
    reg_c = corner(rho,eta);
    rho_c = rho(reg_c);
    eta_c = eta(reg_c);
    return
  end

  % Othersise use local smoothing followed by fitting a 2-D spline curve
  % to the smoothed discrete L-curve. Restrict the analysis of the L-curve
  % according to s_thr.
  if (nargin > 3)
    if (nargin==8)       % In case the bound M is in action.
      s = s(index,:);
    end
    index = find(s > s_thr);
    rho = rho(index); eta = eta(index); reg_param = reg_param(index);
  end

  % Convert to logarithms.
  lr = length(rho);
  lrho = log(rho); leta = log(eta); slrho = lrho; sleta = leta;

  % For all interior points k = q+1:length(rho)-q-1 on the discrete
  % L-curve, perform local smoothing with a polynomial of degree deg
  % to the points k-q:k+q.
  v = (-q:q)'; A = zeros(2*q+1,deg+1); A(:,1) = ones(length(v),1);
  for j = 2:deg+1, A(:,j) = A(:,j-1).*v; end
  for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
  end

  % Fit a 2-D spline curve to the smoothed discrete L-curve.
  sp = spmak((1:lr+order),[slrho';sleta']);
  pp = ppbrk(sp2pp(sp),[4,lr+1]);

  % Extract abscissa and ordinate splines and differentiate them.
  % Compute as many function values as default in spleval.
  P     = spleval(pp);  dpp   = fnder(pp);
  D     = spleval(dpp); ddpp  = fnder(pp,2);
  DD    = spleval(ddpp);
  ppx   = P(1,:);       ppy   = P(2,:);
  dppx  = D(1,:);       dppy  = D(2,:);
  ddppx = DD(1,:);      ddppy = DD(2,:);

  % Compute the corner of the discretized .spline curve via max. curvature.
  % No need to refine this corner, since the final regularization
  % parameter is discrete anyway.
  % Define curvature = 0 where both dppx and dppy are zero.
  k1    = dppx.*ddppy - ddppx.*dppy;
  k2    = (dppx.^2 + dppy.^2).^(1.5);
  I_nz  = find(k2 ~= 0);
  kappa = zeros(1,length(dppx));
  kappa(I_nz) = -k1(I_nz)./k2(I_nz);
  [kmax,ikmax] = max(kappa);
  x_corner = ppx(ikmax); y_corner = ppy(ikmax);

  % Locate the point on the discrete L-curve which is closest to the
  % corner of the spline curve.  Prefer a point below and to the
  % left of the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  if (kmax < 0)
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    index = find(lrho < x_corner & leta < y_corner);
    if ~isempty(index)
      [dummy,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
      rpi = index(rpi);
    else
      [dummy,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); rho_c = rho(rpi); eta_c = eta(rpi);
  end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  % The L-curve is differentiable; computation of curvature in
  % log-log scale is easy.

  % Compute g = - curvature of L-curve.
  g = lcfun(reg_param,s,beta,xi,1);

  % Locate the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  [gmin,gi] = min(g);
  reg_c = fminbnd(GGG,...
    reg_param(min(gi+1,length(g))),reg_param(max(gi-1,1)),...
    optimset('Display','off'),s,beta,xi,1); % Minimizer.
  kappa_max = - lcfun(reg_c,s,beta,xi,1); % Maximum curvature.

  if (kappa_max < 0)
    lr = length(rho);
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
  else
    f = s./(s + reg_c);
    eta_c = norm(f.*xi);
    rho_c = norm((1-f).*beta);
    if (m>n), rho_c = sqrt(rho_c^2 + norm(b0)^2); end
  end

else
  error('Illegal method')
end

function g = lcfun(lambda,s,beta,xi,fifth)

% Auxiliary routine for l_corner; computes the NEGATIVE of the curvature.
% Note: lambda may be a vector.  PCH, IMM, Jan. 4, 2008.

% Initialization.
phi = zeros(size(lambda)); dphi = phi; psi = phi; dpsi = phi;
eta = phi; rho = phi;

% Compute some intermediate quantities.
for i = 1:length(lambda)
  if (nargin==4)
    f  = (s.^2)./(s.^2 + lambda(i)^2);
  else
    f  = s./(s + lambda(i));
  end
  cf = 1 - f;
  eta(i) = norm(f.*xi);
  rho(i) = norm(cf.*beta);
  f1 = -2*f.*cf/lambda(i);
  f2 = -f1.*(3-4*f)/lambda(i);
  phi(i)  = sum(f.*f1.*abs(xi).^2);
  psi(i)  = sum(cf.*f1.*abs(beta).^2);
  dphi(i) = sum((f1.^2 + f.*f2).*abs(xi).^2);
  dpsi(i) = sum((-f1.^2 + cf.*f2).*abs(beta).^2);
end

% Now compute the first and second derivatives of eta and rho
% with respect to lambda;
deta  =  phi./eta;
drho  = -psi./rho;
ddeta =  dphi./eta - deta.*(deta./eta);
ddrho = -dpsi./rho - drho.*(drho./rho);

% Convert to derivatives of log(eta) and log(rho).
dlogeta  = deta./eta;
dlogrho  = drho./rho;
ddlogeta = ddeta./eta - (dlogeta).^2;
ddlogrho = ddrho./rho - (dlogrho).^2;

% Let g = curvature.
g = - (dlogrho.*ddlogeta - ddlogrho.*dlogeta)./...
      (dlogrho.^2 + dlogeta.^2).^(1.5);
end
end

function plot_lc(rho,eta,marker,ps,reg_param)
%PLOT_LC Plot the L-curve.
%
% plot_lc(rho,eta,marker,ps,reg_param)
%
% Plots the L-shaped curve of the solution norm
%    eta = || x ||      if   ps = 1
%    eta = || L x ||    if   ps = 2
% as a function of the residual norm rho = || A x - b ||.  If ps is
% not specified, the value ps = 1 is assumed.
%
% The text string marker is used as marker.  If marker is not
% specified, the marker '-' is used.
%
% If a fifth argument reg_param is present, holding the regularization
% parameters corresponding to rho and eta, then some points on the
% L-curve are identified by their corresponding parameter.

% Per Christian Hansen, IMM, 12/29/97.

% Set defaults.
if (nargin==2), marker = '-'; end  % Default marker.
if (nargin < 4), ps = 1; end       % Std. form is default.
np = 10;                           % Number of identified points.

% Initialization.
if (ps < 1 | ps > 2), error('Illegal value of ps'), end
n = length(rho); ni = round(n/np);

% Make plot.
loglog(rho(2:end-1),eta(2:end-1)), ax = axis;
if (max(eta)/min(eta) > 10 | max(rho)/min(rho) > 10)
  if (nargin < 5)
    loglog(rho,eta,marker), axis(ax)
  else
    loglog(rho,eta,marker,rho(ni:ni:n),eta(ni:ni:n),'x'), axis(ax)
    HoldState = ishold; hold on;
    for k = ni:ni:n
      text(rho(k),eta(k),num2str(reg_param(k)));
    end
    if (~HoldState), hold off; end
  end
else
  if (nargin < 5)
    plot(rho,eta,marker), axis(ax)
  else
    plot(rho,eta,marker,rho(ni:ni:n),eta(ni:ni:n),'x'), axis(ax)
    HoldState = ishold; hold on;
    for k = ni:ni:n
      text(rho(k),eta(k),num2str(reg_param(k)));
    end
    if (~HoldState), hold off; end
  end
end
xlabel('residual norm || A x - b ||_2')
if (ps==1)
  ylabel('solution norm || x ||_2')
else
  ylabel('solution semi-norm || L x ||_2')
end
title('L-curve')
end
end
end