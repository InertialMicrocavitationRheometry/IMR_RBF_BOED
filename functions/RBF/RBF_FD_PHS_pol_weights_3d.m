function w = RBF_FD_PHS_pol_weights_3d (x,y,z,xe,ye,ze,m,d)

% based on codes provided by Flyer et al. (2016)

% Input parameters
% x,y Column vectors with stencil node locations; approximation to
% be accurate at xe,ye
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
%
% Output parameter
% w Matrix with seven columns, containing weights for d/dx, d/dy, d/dz,
% and the Laplacian d2/dx2+d2/dy2
% and d2/dx2, d2/dy2, d2/dz2, 
% , respectively.

x = x-xe; y = y-ye;   z = z-ze; % Shift nodes so stencil centered at origin

n = length(x);

% scale the distance to O(1)

scale = sqrt(x(end)^2+y(end)^2 +z(end)^2 );

x = x/scale; y = y/scale;   z = z/scale;


r =  hypot(hypot(x,y),z);

% ------ RBF part --------------------------------------------------------
A0 = hypot(hypot(bsxfun(@minus,x,x'),bsxfun(@minus,y,y')), bsxfun(@minus,z,z')).^m; % RBF matrix

% L0 = m*(bsxfun(@times,(hypot(x,y)).^(m-2),[-x,-y,m*ones(n,1)])); % RHSs




L0 = m*(bsxfun(@times, r.^(m-2),[-x,-y,-z , (m+1)*ones(n,1), 1+(m-2)*x.^2./(r).^(2),1+(m-2)*y.^2./(r).^(2),1+(m-2)*z.^2./(r).^(2) ])); % RHSs

if x(1) == 0 && y(1) ==0  && z(1) ==0
    L0(1,4)=0;
    L0(1,5) =0;
    L0(1,6) =0;
     L0(1,7) =0;

end

% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
    A = A0; L = L0; % i.e. pure RBF

else % Create matrix with polynomial terms and matching constraints
    X = x(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = y(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    Z = z(:,ones(1,d+1)); Z(:,1) = 1; Z = cumprod( Z,2);

    np = (d+1)*(d+2)*(d+3)/6; % Number of polynomial terms
    XY = zeros(n,np);

    col = 1; % Assemble polynomial matrix block

    for k = 0:d

        for k1 = k:-1:0

            %             XY(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
            for k2 =    k-k1:-1:0

                k3 = k-k1-k2;

                XY(:,col) = X(:,k1+1).*Y(:,k2+1).*Z(:,k3+1);

                col = col+1;
            end

        end

    end


    L1 = zeros(np,7); % Create matching RHSs

    if d >= 1; L1(2,1) = 1; L1(3,2) = 1; L1(4,3) = 1; end

    if d >= 2; L1(5,4) = 2; L1(5,5)=2; L1(8,4) = 2; L1(8,6)=2;  L1(10,4) = 2; L1(10,7)=2;  end

    A = [A0,XY;XY',zeros(col-1)]; % Assemble linear system to be solved

    L = [L0;L1]; % Assemble RHSs
end
% ------ Solve for weights -----------------------------------------------

        W = A\L;

w = W(1:n,:);% Extract the RBF-FD weights

w(:,1:3)= w(:,1:3)/scale;
w(:,4:7)= w(:,4:7)/scale^2;

