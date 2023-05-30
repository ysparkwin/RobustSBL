function [ Ilocs, report ] = SBL_v5p12( A , Y , method , upar , options, Npeaks, peak_sep_in_bins )
% function [ gamma , report ] = SBL_v4( A , Y , Nfreq, Nsource, flag )
% The idea behind SBL is to find a diagonal replica 'covariance' Gamma.
% Minimizing (YY^T / AGA^T + penality) should lead to the correct
% replica selection (up to a bogus scale factor/amplitude).
%
% Inputs
%
% A - Multiple frequency augmented dictionary < n , m , f >
%     n: number of sensors
%     m: number of replicas
%     f: number of frequencies
%   Note: if f==1, A = < n , m >
% Y - Multiple snapshot multiple frequency observations < n , L , f >
%     n: number of sensors
%     L: number of snapshots
%     f: number of frequencies
%   Note: if f==1, Y = < n , L >
% options - see SBLset.m 
%
%
% Outputs
%
% gamma <m , 1> - vector containing source power
%                 1: surfaces found by minimum error norm
% report - various report options
%
%--------------------------------------------------------------------------
% Version 1.0:
% Code originally written by P. Gerstoft.
%
% Version 2.23
% Edited to include multiple frequency support: 5/16/16
%
% Version 3.1
% Different convergance norm and code update
% A and Y have now one more dimensions
% Posterior unbiased mean
% handles single snapshot
%
% Version 3.32
% Robust version of SBL with Phi parameter in options structure
%
% Santosh Nannuru & Kay L Gemba
% NoiseLab/SIO/UCSD gemba@ucsd.edu & snannuru@ucsd.edu

%%
options.SBL_v = '5.12';
activeIndRepN = options.activeIndRepN;

%% slicing
if options.tic == 1
    tic
end

%% Initialize variables
Nfreq       = size(A,3);        % number of frequencies
Nsensor     = size(A,1);        % number of sensors
Ntheta      = size(A,2);        % number of dictionary entries
Nsnapshot   = size(Y,2);        % number of snapshots in the data covariance
Nsource     = options.Nsource;  % number of sources

robustSBL = true; % default is to use robust SBL-T method 
if strcmpi(method,'SBL-G')
    options.SBL_v = [options.SBL_v,'-G'];
    ufun = @(t,v) ones(size(t)); % this is actually not used 
    b = 1;
    upar = 1;
    robustSBL = false;
elseif (strcmpi(method,'SBL-T') || isempty(method) )
    options.SBL_v = [options.SBL_v,'-T'];
    ufun = @(t,v) (v+2*Nsensor)./(v + 2*t);  % weight function is t-loss 
    if nargin < 5
       upar = 2.1; % d.o.f v=2.1 is used as the default parameter for t-loss
    end   
    assert( isreal(upar) && upar>=2 && ~isinf(upar),'loss parameter is incorrect');  
    b = tloss_consistency_factor(Nsensor,upar,false); % concistency factor    
elseif (strcmpi(method,'SBL-H') || isempty(method) )   
    options.SBL_v = [options.SBL_v,'-H'];
    ufun = @(t,c) ((t<=c) + (c./t).*(t>c)); % weight function u(t)
    if nargin < 5
       q = 0.85;
    else
       q = upar;
    end          
    assert(isreal(q) && isfinite(q) && (0<q) && (q <1),'loss parmeter is incorrect'); 
    upar = chi2inv(q,2*Nsensor)/2;  
    b = chi2cdf(2*upar,2*(Nsensor+1))+(upar/Nsensor)*(1-q);    
elseif  strcmpi(method,'SBL-Tyl')
    options.SBL_v = [options.SBL_v,'-Tyl'];
    ufun = @(t,c) c./t; 
    [~,wei] = TylM(Y.');
    if nargin < 5
       upar = Nsensor;
    end
    b = 1/median(1./wei);
else
    error('please specify method as a string equalto SLB-G, SBL-T, SBL-H, or SBL-Tyl\n');
end
const = 1/(b*Nsnapshot);

% posterior
if isfield(options,'x_post_calculation') == 1
x_post    = zeros( Ntheta, Nsnapshot,Nfreq);
end
% minimum (global) gamma
gmin_global = realmax;
% space allocation
errornorm   = zeros(options.convergence.maxiter,1);
RYerrornorm = zeros(options.convergence.maxiter,1);
if options.activeIndices == 1
activeInd   = zeros(options.convergence.maxiter,Nsource);
sigcItr   = zeros(options.convergence.maxiter,1);
end

% initialize with CBF output (assume single frequency)
gamma        = zeros(Ntheta, 1);
sigc         = zeros(Nfreq , 1);
for iF = 1:Nfreq
    Af = squeeze(A(:,:,iF));
    Yf = squeeze(Y(:,:,iF));
    Ryyf = (1/Nsnapshot) * (Yf * Yf');

    gamma_init = sum(abs(Af' * Yf).^2, 2) / Nsnapshot ./ sum(abs(Af(:,1)).^2)^2; % Note: normalize with norm(Af(:,iColumn))^4, Below Eq. (24) [Initialization] --*

    % -- Compute the initial noise variance
    if Nsource==1
        [~,indx] = max(gamma_init);
    else
        [~,indx] = SBLpeaks_1D(gamma_init,Nsource);
    end
    Am = Af(:,indx);        % only active replicas
    Aplus = pinv(Am);
    P_N = eye(Nsensor)-Am*Aplus;
    sigc(iF) = real(trace(P_N*Ryyf))/(Nsensor - Nsource);  % noise covariance estimate

    % -- Compute initial estimates of gamma
    gamma_temp =  max(  max(gamma_init)*options.gamma_range , gamma_init -  sigc(iF)); % Note: norm(Af(:,iColumn)) = sqrt( Nsensor ) Below Eq. (24) [Initialization] --*
    gamma = gamma + gamma_temp;
end

% Sample Covariance Matrix
SCM      = zeros( Nsensor , Nsensor , Nfreq );
maxnoise = zeros( Nfreq , 1 );
for iF = 1 : Nfreq
    SCM(:,:,iF) = squeeze(Y(:,:,iF)) * squeeze(Y(:,:,iF))' / Nsnapshot;
    maxnoise(iF)= real(trace(squeeze(SCM(:,:,iF))))/Nsensor;
end

RY           = zeros( Nsensor , Nsensor , Nfreq );
gamma_num    = zeros(Ntheta,Nfreq);
gamma_denum  = zeros(Ntheta,Nfreq);

%% Main Loop
disp(['SBL version ', options.SBL_v ,' initialized.']);

for j1 = 1 : options.convergence.maxiter
        
    % for error analysis
    gammaOld    = gamma;
    RYOld       = RY;
    Itheta      = find(gamma>max(gamma)*options.gamma_range);
    gammaSmall  = gamma(Itheta);

    %% gamma update
    % this is a sum over frequency and can be done by multiple processors
    % --> multi-frequency SBL should be almost as fast as single frequency if num(proc) = Nfreq
    for iF = 1 : Nfreq
        Af = squeeze(A(:,Itheta,iF));
        Yf = squeeze(Y(:,:,iF));

        % SigmaY inverse
        SigmaYinv = eye(Nsensor) / (sigc(iF) * eye(Nsensor) + ...
           Af * bsxfun(@times, gammaSmall, Af') );

        if robustSBL
            t = real(sum(conj(Yf).*(SigmaYinv*Yf))); % norms
            u = ufun(t,upar);
            if strcmpi(method,'SBL-Tyl')
                RY(:,:,iF) = (1/Nsnapshot)*(Yf.*repmat(u,Nsensor,1))*Yf';
%                 RY(:,:,iF) = (1/b)*(RY(:,:,iF)/(real(trace(RY(:,:,iF) ))));
                RY(:,:,iF) = (1/b)*(RY(:,:,iF)/(real(trace(RY(:,:,iF)/Nsensor)))); % Peter's suggestion
            else
                RY(:,:,iF) = const*(Yf.*repmat(u,Nsensor,1))*Yf'; %RY(:,:,iF) = (1/b/Nsnapshot)*(Yf*diag(u))*Yf';
            end
        else, RY(:,:,iF) = SCM(:,:,iF);
        end
        B =  SigmaYinv*Af; % \Sigma^-1 a_m , m=1,..,M

        % Sum over snapshots and normalize
        gamma_num(Itheta,iF)   = ( real( sum ( conj(B).*(RY(:,:,iF)*B)) ) );

        % positive def quantity, abs takes care of roundoff errors        
        gamma_denum(Itheta,iF) = abs( sum  ( (Af' * SigmaYinv).' .* Af, 1 ) )  ;
    end
    
    % Fixed point Eq. update
    gamma(Itheta)  = gamma(Itheta)   .* ((sum( gamma_num(Itheta,:)    ,2 ) ./...
        sum( gamma_denum(Itheta,:)  ,2 ) ).^(1/options.fixedpoint) ) ;

    
    %% sigma update
    % locate same peaks for all frequencies
    if j1 == 1, IlocsOld = [];
    else,       IlocsOld = Ilocs;
    end

    [~, Ilocs] = SBLpeaks_1D(gamma, Nsource);
    sigmaUp = 1;
    if ~robustSBL && isempty(setdiff(IlocsOld,Ilocs)), sigmaUp = 0; end %SBL-G skips sigma update when active indices did not change.    
    if sigmaUp
        for iF = 1 : Nfreq
            % only active replicas
            Am       = squeeze(A(:,Ilocs,iF));
            % noise estimate
            sigc(iF) = real(trace( (eye(Nsensor)-Am*pinv(Am)) * squeeze(RY(:,:,iF)) ) / ( Nsensor - Nsource ) );
            sigc(iF) = min(sigc(iF),maxnoise(iF));          % cannot be larger than signal+noise.
            sigc(iF) = max(sigc(iF),maxnoise(iF)*10^-10);   % snr>100 is unlikely larger than signal.
        end
    end

    % check if the index repeats
    if ~isempty(setdiff(IlocsOld,Ilocs))
        activeIndRep = 0;
    else
        if exist('activeIndRep','var')==0, activeIndRep = 0; end
        activeIndRep = activeIndRep + 1;
    end

    if options.activeIndices == 1
    activeInd(j1,:) = Ilocs.';
    sigcItr(j1)     = sigc;
    end    
        
    %% Convergence checks convergance and displays status reports
    % convergence indicator
    errornorm(j1)   = norm ( gamma - gammaOld, 1 ) / norm ( gamma, 1 );
    RYerrornorm(j1) = norm ( RY - RYOld, 'fro' ) / norm ( RY ,'fro' );
    report.results.gamma(j1)= max(gamma);
    % global min error
    if j1 > options.convergence.min_iteration  &&  errornorm(j1) < gmin_global
        gmin_global  = errornorm(j1);
        gamma_min    = gamma;
        iteration_L1 = j1;
    end
    
    % inline convergence code
    if activeIndRep > activeIndRepN || (j1 > options.convergence.min_iteration && ( errornorm(j1) < options.convergence.error  || iteration_L1 + options.convergence.delay <= j1))
        if options.flag == 1
            disp(['Solution converged. Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
        break; % goodbye
    % not convereged
    elseif j1 == options.convergence.maxiter
        if options.flag == 1
            warning(['Solution not converged. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
    % status report
    elseif j1 ~= options.convergence.maxiter  && options.flag == 1 && mod(j1,options.status_report) == 0 % Iteration reporting
        disp(['Iteration: ',num2str(sprintf('%.4u',j1)),' Dic size: ',num2str(length(Itheta)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.' ])
    end
    
end

%% Peak finding
[~, Ilocs] = SBLpeaks_1D_sep(gamma, Npeaks, peak_sep_in_bins);

%% Posterior distribution for polarity
% x_post - posterior unbiased mean
if isfield(options,'x_post_calculation') == 1
for iF = 1 : Nfreq
    Af = squeeze (A(:,:,iF));
    x_post(:,:,iF) = repmat(gamma, [1 Nsnapshot] ) .* (Af' / (sigc(iF) * eye(Nsensor) + ...
        Af * (repmat(gamma, [1 Nsensor] ) .* Af')) * squeeze(Y(:,:,iF)));
end
end
%% function return
gamma = gamma_min;

%% Report section
% vectors containing errors
report.results.error    = errornorm(1:iteration_L1);
report.results.RYerror  = RYerrornorm(1:iteration_L1);
% Error when minimum was obtained
report.results.iteration_L1 = iteration_L1;
% General info
report.results.final_iteration.iteration  = j1;
report.results.final_iteration.noisepower = sigc;

if options.tic == 1
    report.results.toc = toc;
else
    report.results.toc = 0;
end

% data
report.results.final_iteration.gamma  = gamma;
report.results.final_iteration.RY  = RY;
if isfield(options,'x_post_calculation') == 1
report.results.final_iteration.x_post = x_post;
end

report.options = options;

if options.activeIndices == 1
report.results.activeIndices = activeInd(1:iteration_L1,:);
report.results.sigcItr       = sigcItr(1:iteration_L1);
end

% debug output parameters (assuming single frequency)
%report.SigmaYinv = SigmaYinv;
%report.SCM = squeeze(SCM);
end

function [pks, locs] = SBLpeaks_1D(gamma, Nsources)
% fast alternative for findpeaks in 1D case
% output variables
pks  = zeros(Nsources,1);
locs = zeros(Nsources,1);

% zero padding on the boundary
gamma_new           = zeros(length(gamma)+2,1);
gamma_new(2:end-1)  = gamma;

[~, Ilocs] = sort(gamma,'descend');

% current number of peaks found
npeaks = 0;
for ii = 1:length(Ilocs)
    % local patch area surrounding the current array entry i.e. (r,c)
    local_patch     = gamma_new(Ilocs(ii):Ilocs(ii)+2);
    % zero the center
    local_patch(2)  = 0;
    if sum(sum(gamma(Ilocs(ii)) > local_patch)) == 3
        npeaks       = npeaks + 1;
        pks(npeaks)  = gamma(Ilocs(ii));
        locs(npeaks) = Ilocs(ii);
        if npeaks == Nsources % if found sufficient peaks, break
            break;
        end
    end
end
if npeaks ~= Nsources % if Nsources not found
    pks(npeaks+1:Nsources)  = [];
    locs(npeaks+1:Nsources) = [];
end
end

function [pks, locs] = SBLpeaks_1D_sep(gamma, Nsources, NIndSep)
% fast alternative for findpeaks in 1D case with an index separation
% output variables
pks  = zeros(Nsources,1);
locs = zeros(Nsources,1);

% zero padding on the boundary
gamma_new           = zeros(length(gamma)+2,1);
gamma_new(2:end-1)  = gamma;

[~, Ilocs]= sort(gamma,'descend');

% current number of peaks found
npeaks = 0;
for K = 1:Nsources
    for ii = 1:length(Ilocs)
        % local patch area surrounding the current array entry i.e. (r,c)
        local_patch     = gamma_new(Ilocs(ii):Ilocs(ii)+2);
        % zero the center
        local_patch(2)  = 0;
        if sum(sum(gamma(Ilocs(ii)) > local_patch)) == 3
            npeaks       = npeaks + 1;
            pks(npeaks)  = gamma(Ilocs(ii));
            locs(npeaks) = Ilocs(ii);
            Ilocs( abs(Ilocs - Ilocs(ii)) <= NIndSep ) = [];
            break;
        end
        if ii == length(Ilocs)
            error('It was not possible to localize peaks with the chosen sep!')
        end
    end
    if isempty(Ilocs)
        error('It was not possible to localize peaks with the chosen sep!')
    end
    if npeaks == Nsources % if found sufficient peaks, break
        break;
    end
end
end

function b=tloss_consistency_factor(p,v,realdata)
%  Function that computes the scaling factor for multivariate t-weight
%  function so that the returned scatter matrix is concistent estimator of
%  the covariance matrix under the assumption that the data is from 
%  Gaussian distribution
%--------------------------------------------------------------------------   

% First try by numerical integration 
if realdata 
     b = tloss_consistency_factor_int(p,v);
else
     b = tloss_consistency_factor_int(2*p,v);        
end

% If integration did not converge, then compute b by MC simulation 
if isnan(b)
   % need to use MC simul to find b
   MCsimul = 100000;
   if realdata
        t = chi2rnd(p,1,MCsimul);
        psifun = @(t,v) (p+v)*t./(v+t); % weight function
   else 
        t = (1/2)*chi2rnd(2*p,1,MCsimul);
        psifun = @(t,v) (v+2*p)*t./(v+2*t); % weight function
   end
   b = (1/p)*mean(psifun(t,v));
end
end

function b=tloss_consistency_factor_int(p,v)
% computes the concistency factor b = (1/p) E[|| x ||^2 u_v( ||x||^2)] when
% x ~ N_p(0,I). 
sfun = @(x,p,v)  (x.^(p/2)./(v+ x) .* exp(-x/2));
c = 2^(p/2)*gamma(p/2);
w = warning('off','all');
q = (1/c)*integral(@(x)sfun(x,p,v),0,Inf);
b = ((v+p)/p)*q; % consistency factor  
warning(w)
end

function [C,wei,invC,iter] = TylM(X,varargin)
% TylM compute the M-estimator of scatter using Huber's weight function and
% the fixed-point (FP) algorithm. 
%  
% Data is assumed to be centered (or the symmetry center parameter = 0)
%   
% USAGE
% C = TylM(X);
% [C,invC] = TylM(X,'invC',cov(X) \ eye(p),'printitn',1);
% 
% INPUT
%   X       : the data matrix with n rows (observations) and p columns.
%             observations can be real or complex-valued. 
%
% Name-Value Pair Arguments
% -------------------------
% TylM can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% Name      Value and description
%==========================================================================
%  Optional inputs: 
% 'invC' :  [] (default) or a positive definite  p x p matrix
%       correspond to the inverse scatter matrix to start the FP iterations.  
%       If not specified ([]), algorithm uses the  inverse of the SCM as the
%       initial start for the FP algorithm.
%   
% 'printitn' : non-negative integer (Default is 0) 
%       Print iteration convergence. If 1, then print each iteration; if 2, 
%       then every 2nd iteration, etc. 
% 
% OUTPUT
%   C       : Tyler's M-estimator of scatter 
%   invC    : the inverse matrix of C
%   iter    : number of iterations of the FP algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

argin = inputParser;
[n, p] = size(X);

valX = @(X) assert(numel(size (X))==2 && size(X,1) >= size(X,2) && size(X,2) >= 2, ... 
        '''X'' must be a matrix having at least two columns (variables) and n >= p');
addRequired(argin,'X',valX);

if any (any (isnan (X)))
    error ('Input data contains NaN''s.');
end

if ~isa (X, 'double')
    X = double(X);
end

realdata = true;
if ~isreal(X)
   realdata=false; 
end

%-- optional input parsing rules
valinvC = @(x) assert( isempty(x) || ( ismatrix(x) && size(x,2)==size(x,1) && size(x,1)==p), ...
    'the input invC must be a empty set or matrix of size p x p');
addParameter(argin,'invC',[],valinvC);

valprintitn = @(x) assert( isscalar(x) && x >= 0 && (x==round(x)),  ... 
    '''printitn'' must be a non-negative integer');
addParameter(argin,'printitn',0, valprintitn);

%---  parse inputs
parse(argin,X,varargin{:});
invC        = argin.Results.invC;
printitn    = argin.Results.printitn; 

if p > n
    error(['data matrix X: sample size (# of rows) exceeds the data ' ...
        ' dimensionality (# of columns)']);
end

if isempty(invC) 
    C = X'*X/n;
    C = p*C/trace(C);
    invC = C\eye(p);
end
           
MAX_ITER = 5000; % Max number of iteration
EPS = 5.0e-4;    % Iteration accuracy
iter = 1;

while (iter<MAX_ITER)
      t = real(sum((X*invC).*conj(X),2)); % norms        
      C = (p/n)*X'*(X.*repmat(1./t,1,p)); 
      C = p*C/trace(C);
      err = eye(p)-invC*C;
      d = norm(err(:),Inf);
      if mod(iter,printitn)==0
            fprintf('At iter = %4d, dis=%.6f\n',iter,d);
      end
      invC = C\eye(p);
      if (d<=EPS), break; end         
      iter = iter+1;  
end
C(1:(p+1):p^2) = real(C(1:(p+1):p^2)); 
wei = p./t;
if (iter==MAX_ITER)
     error(['WARNING! Slow convergence: error of the solution is %f\n ' ...
         'after %d iterations\n'],d,iter);
end
end