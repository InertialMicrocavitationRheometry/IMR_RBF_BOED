
function w = RBF_interp_PHS_pol_weights(X0,XE0,m,d)

% based on codes provided by Flyer et al. (2016)

% Input parameters
% x,y Column vectors with stencil node locations; approximation to
% be accurate at xe,ye
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials
%
% Output parameter
% w Matrix with three columns, containing weights for 1.


x = X0(:,1); y = X0(:,2); 
xe = XE0(:,1); ye = XE0(:,2);

x = x-xe; y = y-ye; % Shift nodes so stencil centered at origin

n = length(x);

% dr  =  mean(sqrt( (x(1)-x(2:5)).^2+(y(1)-y(2:5)).^2));

dx  =  mean(sqrt( (x(1)-x(2:5)).^2));
dy  =  mean(sqrt( (y(1)-y(2:5)).^2));
% 
 dr = sqrt(dx^2+dy^2);

% scale the distance to O(1)

scale = sqrt(x(end)^2+y(end)^2);

  % x = x/scale; y = y/scale;

A0 = hypot(bsxfun(@minus,x,x'),bsxfun(@minus,y,y')).^m; % RBF matrix

L0 = (bsxfun(@times,(hypot(x,y)).^(m),1)); % RHSs



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

     L1 = zeros(np,1); % Create matching RHSs

     L1(1) =1;
    A = [A0,XY;XY',zeros(col-1)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs


end
% ------ Solve for weights -----------------------------------------------

            W = A\L;

    % [Yu,Yd,Yv]          =    svd(A);
    % epsilon             =    1e-6;
    % Yd                  =    diag(Yd);
    % [Yd,sort_i]         =    sort(Yd,'descend');
    % Yv                  =    Yv(:,sort_i);
    % Yv                  =    Yv(:,Yd>epsilon);
    % Yu                  =    Yu(:,sort_i);
    % Yu                  =    Yu(:,Yd>epsilon);
    % Yd                  =    Yd(Yd>epsilon);
    % 
    % Ainv                =    Yv*diag(1./Yd)*Yu';
    % W                   =    Ainv*L;


w = W(1:n,:);% Extract the RBF weights

  % w(:,1)= w(:,1)/scale;

