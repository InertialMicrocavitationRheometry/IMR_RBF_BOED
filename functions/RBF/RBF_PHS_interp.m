function [D0] = RBF_PHS_interp(xy1,xy_s,Nearest_Idx_interp,k,opt,m,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input

%   xy1:                    approximate derivatives at xy1
%   xy_s:                   given function values at xy_s
%   Nearest_Idx_interp:     Indices of nearest nodes
%   k:                      stencil size
%   opt:                    linear operator, e.g. {1},{x},{y},{L}
%   m:                      order of PHS-RBFs
%   d:                      degree of polynomials

%% Output

%  D0:                      Global differentiation matrix from xy_s to xy1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% References: 
% [1] Chu, T., & Schmidt, O. T. (2023). RBF-FD discretization of the Navier-Stokes 
% equations on scattered but staggered nodes. 
% Journal of Computational Physics, 474, 111756.
% [2] Chu, T., & Schmidt, O. T. (2024). Mesh-free hydrodynamic stability. 
% Journal of Computational Physics, 502, 112822.

%%


lx = max(xy_s(:,1))- min(xy_s(:,1));

ly = max(xy_s(:,2))- min(xy_s(:,2));


if lx<1e-6
    lx=1e-6;
end

 if ly<1e-6
    ly=1e-6;   
 end

xy_s(:,1) = xy_s(:,1)/lx;
xy_s(:,2) = xy_s(:,2)/ly;

xy1(:,1) = xy1(:,1)/lx;
xy1(:,2) = xy1(:,2)/ly;


%%

X1            =           xy1(:,1);
NxNy1         =           length(X1);
X2            =           xy_s(:,1);
NxNy2         =           length(X2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m1 = 1:NxNy1

    if mod(m1,100)==0
        disp(['j= ' num2str(m1)])
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k0 =k;

    xk1 = xy_s(Nearest_Idx_interp(m1,1:k0),:);     %  local stencil

    xe1 = xy1(m1,:);                               %  evaluation point

    % move center to the origin


    xk1 = xk1-xe1;
    xe0 = xe1;
    xe1 = xe1-xe1;


    if opt == '1'

        w = RBF_interp_PHS_pol_weights (xk1,xe1,m,d);      % RBF-PHS method

        weight(m1,1:k0)  =  w(:,1);

    else
        w = RBF_FD_PHS_pol_weights (xk1,xe1,m,d);      % RBF-PHS method

        if  length(opt)==1 && opt =='x'
            weight(m1,1:k0)  =  w(:,1)/lx;
        elseif length(opt)==1 && opt == 'y'
            weight(m1,1:k0)  =  w(:,2)/ly;
        elseif opt == 'L'
            weight(m1,1:k0)  =  w(:,3)/(lx^2+ly^2);
        elseif length(opt)==2 && opt(1) == 'x' &&  opt(2) == 'x'
            weight(m1,1:k0)  =  w(:,4)/lx^2;
        elseif  length(opt)==2 && opt(1) == 'y' &&  opt(2) == 'y'
            weight(m1,1:k0)  =  w(:,5)/ly^2;
        elseif  length(opt)==2 && opt(1) == 'x' &&  opt(2) == 'y'
            weight(m1,1:k0)  =  w(:,6)/lx/ly;
        elseif opt == 'n'
            weight(m1,1:k0)  =  (xe0(1)*w(:,1)+xe0(2)*w(:,2))/sqrt(xe0(1)^2+xe0(2)^2);
        end

        %
    end

end


W                                    =       reshape(weight',1,[]);

Idx_x                                  =       reshape((1:NxNy1).*ones(k,1),1,[]);
Idx_y                                  =       reshape(Nearest_Idx_interp(:,1:k)',1,[]);


D0                               =       sparse(Idx_x,Idx_y,W,NxNy1,NxNy2);


end