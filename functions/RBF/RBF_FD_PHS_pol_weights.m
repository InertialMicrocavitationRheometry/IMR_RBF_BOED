

function [w] = RBF_FD_PHS_pol_weights(X0,XE0,m,d,varargin)

% based on codes provided by Flyer et al. (2016)

% Input parameters
% x,y Column vectors with stencil node locations; approximation to be accurate at xe,ye
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
%
% Output parameter
% w Matrix with six columns, containing weights for d/dx, d/dy,the Laplacian d2/dx2+d2/dy2,d2/dx2, d2/dy2, and d2/dxy, respectively.


x = X0(:,1); y = X0(:,2); 
xe = XE0(:,1); ye = XE0(:,2);

if nargin==5

    scale = varargin{1};
else
    scale = sqrt(x(end)^2+y(end)^2)/1;
end

scale_x = scale; scale_y = scale;

%  scale_x = max(x)-min(x);
%  
%  scale_y = max(y)-min(y);

x = x-xe; y = y-ye; % Shift nodes so stencil centered at origin
n = length(x);

 x = x/scale_x; y = y/scale_y;



% ------ RBF part --------------------------------------------------------
A0 = hypot(bsxfun(@minus,x,x'),bsxfun(@minus,y,y')).^m; % RBF matrix

% L0 = m*(bsxfun(@times,(hypot(x,y)).^(m-2),[-x,-y,m*ones(n,1)])); % RHSs


L0 = m*(bsxfun(@times,(hypot(x,y)).^(m-2),[-x,-y,m*ones(n,1),1+(m-2)*x.^2./(hypot(x,y)).^(2),1+(m-2)*y.^2./(hypot(x,y)).^(2),(m-2)*x.*y./(hypot(x,y)).^(2)])); % RHSs


for j = 1:n

    if x(j) == 0 && y(j) == 0
        L0(j,4) =0;
        L0(j,5) =0;
        L0(j,6) =0;
    end

end

% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
    A = A0; L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
    X = x(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = y(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    np = (d+1)*(d+2)/2; % Number of polynomial terms
    XY = zeros(n,np); col = 1; % Assemble polynomial matrix block
    for k = 0:d
        XY(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
        col = col+k+1;
    end
    L1 = zeros(np,6); % Create matching RHSs
    if d >= 1; L1(2,1) = 1; L1(3,2) = 1; end
    if d >= 2; L1(4,3) = 2; L1(4,4)=2; L1(6,3) = 2; L1(6,5)=2; L1(5,6)=1;end
    A = [A0,XY;XY',zeros(col-1)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs
end
% ------ Solve for weights -----------------------------------------------
%  

            A = (A+A')/2;
% 
               W = A\L;

     % [U,S,V]   =   svd(A);
     % 
     %        S  =  diag(S);
     % 
     %        id =  (S<1e-8);
     % 
     %        U(:,id)=[]; V(:,id)=[]; S(id)=[];
     %        S      =    diag(1./S);
     %        A_inv =    V*S*(U');
     % 
     %    W = A_inv*L;


w = W(1:n,:);      % Extract the RBF-FD weights


w(:,1)= w(:,1)/scale_x;
w(:,2)= w(:,2)/scale_y;

w(:,3)= w(:,3)/(scale^2);

w(:,4)= w(:,4)/scale_x^2;
w(:,5)= w(:,5)/scale_y^2;
w(:,6)= w(:,6)/scale_x/scale_y;
