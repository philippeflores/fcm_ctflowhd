function [y,lambda,erreur3,erreur1,flag,R] = PCTF3D(dataMarg,R,varargin)
% AOADMMCOUPLE (Alternate Optimization using Alternating Direction Method of 
% Multipliers for Coupled Tensor Factorization) with Sum To One constraints
% 
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%
% dataX : structure containing data
% R : Rank of factorization
% triplets : indices of each triplets considered for coupled TF
% N1 : Maximal number of outer iterations
% N2 : Maximal number of inner iterations
% eps : tolerance for stopping criteria of the inner loop
% y0 : initialisation
% 
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% 
% y : 1xM cell that contains the factors of the factorization of X (where N
% is the order of the tensor)
% lambda : Rx1 vector containing the proportion of each components on the
% decomposition
% erreur : Vector that contains frobenius errors between the factorization
% tensor and the tensor X
% flag : information about the stop of the algorithm
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Set the default parameters
% -------------------------------------------------------------------------
T1 = 1000;
T2 = 20;
eps = 10^(-10);
boolComputeError = 0;
tolNewY = 10^-7;
boolPlot = 1;
% -------------------------------------------------------------------------
% Read the input arguments
% -------------------------------------------------------------------------

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'Y0'
                y0 = varargin{i+1};
			case 'LAMBDA0'
				lambda0 = varargin{i+1};
            case 'T1'
                T1 = varargin{i+1};
            case 'T2'
                T2 = varargin{i+1};
            case 'EPS'
                eps = varargin{i+1};
            case 'TOLNEWY'
                tolNewY = varargin{i+1};
            case 'BOOLPLOT'
                boolPlot = varargin{i+1};
            case 'THRESHLAMBDA'
                threshLambda = varargin{i+1};
            case 'COMPUTEERROR'
                boolComputeError = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

rho = -1;
M = dataMarg.M;
I = dataMarg.I;
calT = dataMarg.calT;
marg = dataMarg.marg;
indVar = dataMarg.indVar;

y = cell(1,M); yt = cell(1,M);

if exist('y0','var')
	for m = 1:M, y{m} = y0{m}; end
else
	for m = 1:M, y{m} = rand(I,R); end
end
for m = 1:M
	yt{m} = zeros(size(y{m}));
	for r = 1:R, y{m}(:,r) = y{m}(:,r)/sum(y{m}(:,r)); end
end

if exist('lambda0','var'), lambda = lambda0; else, lambda = rand(R,1); end
lambda = lambda/sum(lambda); lambdat = zeros(R,1);

erreur3 = inf(T1+1,1);
erreur1 = inf(T1+1,1);
[erreur3(1),erreur1(1)] = computeError(y,lambda,marg,calT,indVar);

flag = 0; t1 = 0;
if boolPlot==1, fprintf("\nPCTF3D :    "), end
while flag == 0
    if boolPlot==1
	    for i = 1:99
		    trigPerc = floor(i*T1/100)+1;
		    if t1==trigPerc, fprintf('\b\b\b\b %2d%%',i); end
	    end
    end

    for m = 1:M

		tempy= y;
		[y{m},yt{m}] = ADMMy(y,yt,m,lambda,marg,calT,R,rho,eps,T2,indVar);
		
    end
    
    tempLambda = lambda;
	[lambda,lambdat] = ADMMl(lambda,lambdat,y,marg,calT,R,rho,eps,T2,indVar);

	if boolComputeError==1
		[erreur3(t1+1),erreur1(t1+1)] = computeError(y,lambda,marg,calT,indVar);
	end
	t1 = t1+1;
	if t1>T1, flag = T1; end
    distNewY = 0;
    for t = 1:size(calT,2)
        posVar = findPosVar(indVar,calT{t});
        distNewY = distNewY+frob(cpdgen({y{posVar(1)},y{posVar(2)},y{posVar(3)}},lambda)-cpdgen({tempy{posVar(1)},tempy{posVar(2)},tempy{posVar(3)}},tempLambda));
    end
    tempOld = [reshape(vertcat(tempy{:}), I*R*M, 1); tempLambda];
    tempNew = [reshape(vertcat(y{:}), I*R*M, 1); lambda];
    if norm(tempOld - tempNew) <= tolNewY
        flag = t1+tolNewY;
    end

end

[lambda,indLambda] = sort(lambda,'descend');
y = permfacto(y,indLambda);

if boolPlot==1, fprintf('\b\b\b\b 100%%'), end

if exist('threshLambda','var')
	for r = R:-1:1, if lambda(r)<threshLambda, for m = 1:M, y{m}(:,r) = []; end, end, end

	if size(y{1},2)<R
        if boolPlot==1
            if R-size(y{1},2)==1
		        fprintf("\nOne component has been deleted. New rank is %d.\n",size(y{1},2));
            else
		        fprintf("\n%d components have been deleted. New rank is %d.\n",R-size(y{1},2),size(y{1},2));
            end
        end
		R = size(y{1},2); lambda = lambda(1:R);
	end
	[erreur3(end),erreur1(end)] = computeError(y,lambda,marg,calT,indVar);
end

erreur3(isinf(erreur3)) = [];
erreur1(isinf(erreur1)) = [];
end

function [H,U] = ADMMy(y,yt,m,lambda,marg,calT,R,rho,eps,T2,indVar)

G = zeros(R); V = zeros(size(yt{m}'));
for nTrip = 1:size(calT,2)
	posVar = findPosVar(indVar,calT{nTrip});
	if any(posVar==m)
		indKrProd = posVar;
		indn = find(posVar==m);
		indKrProd(indn) = [];
		krProd = krb(y{indKrProd(2)},y{indKrProd(1)});
		G = G+krProd'*krProd;
		V = V+krProd'*nshape(marg{nTrip},[indn 1:indn-1 indn+1:3])';
	end
end

G = G.*(lambda*lambda');
V = diag(lambda)*V;

if rho<0, rho = trace(G/R); end
L = pinv(G+rho*eye(R));

innerFlag = 0;

t2 = 1;

H = y{m};
U = yt{m};
Ht = zeros(size(H'));

while innerFlag ==0
	
	Htm = Ht;
	Ht = L*(V+rho*(H+U)');
	
	H = ProjectOntoSimplex(Ht'-U,1);
		
	U = U+H-Ht';
	
	r = frob(H-Ht')^2;
	s = frob(rho*(Ht-Htm))^2;
	
	% Incrementation of the count
	t2 = t2+1;
	
	% Update of the inner flag
	if t2>T2, innerFlag = T2; end
	if r<eps && s<eps, innerFlag = eps; end
		
end

end

function [h,u] = ADMMl(lambda,lambdat,y,marg,calT,R,rho,eps,T2,indVar)

% Computation of G
G = zeros(R); V = zeros(size(lambda));
for nTrip = 1:size(calT,2)
	posVar = findPosVar(indVar,calT{nTrip});
	A = y{posVar(1)}; B = y{posVar(2)}; C = y{posVar(3)};
	G = G+(A'*A).*(B'*B).*(C'*C);
	V = V+krb(C,krb(B,A))'*marg{nTrip}(:);
end

if rho<0, rho = trace(G)/R; end
L = pinv(G+rho*eye(R));

innerFlag = 0;

t2 = 1;

h = lambda;
u = lambdat;
ht = zeros(size(lambda));

while innerFlag ==0
	
	htm = ht;
	ht = L*(V+rho*(h+u));
	
	h = ProjectOntoSimplex(ht-u,1);
	
	u = u+h-ht;
	
	r = frob(h-ht')^2;
	s = frob(rho*(ht-htm))^2;
	
	t2 = t2+1;
	
	if t2>T2, innerFlag = T2; end
	if r<eps && s<eps, innerFlag = eps; end
	
end

end
