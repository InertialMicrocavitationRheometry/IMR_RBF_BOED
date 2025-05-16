function [D0] = RBF_PHS_interp_4d(xy1,xy_s,Nearest_Idx_interp,k,opt,m,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input

%   xy1:                    approximate derivatives at xy1
%   xy_s:                   given function values at xy_s
%   Nearest_Idx_interp:     nearest nodes
%   k:                      stencil size
%   opt:                    linear operator, e.g. {x},{y},{L}
%   m:                      order of PHS-RBFs
%   d:                      degree of polynomials

%% Output

%  D0:                      Global differentiation matrix from xy_s to xy1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lx = max(xy_s(:,1))- min(xy_s(:,1));
ly = max(xy_s(:,2))- min(xy_s(:,2));
lz = max(xy_s(:,3))- min(xy_s(:,3));
lx4 = max(xy_s(:,4))- min(xy_s(:,4));


if lx<1e-6
    lx=1e-6;
end

 if ly<1e-6
    ly=1e-6;   
 end

  if lz<1e-6
    lz=1e-6;   
  end

  if lx4<1e-6
    lx4=1e-6;
end

xy_s(:,1) = xy_s(:,1)/lx;
xy_s(:,2) = xy_s(:,2)/ly;
xy_s(:,3) = xy_s(:,3)/lz;
xy_s(:,4) = xy_s(:,4)/lx4;


xy1(:,1) = xy1(:,1)/lx;
xy1(:,2) = xy1(:,2)/ly;
xy1(:,3) = xy1(:,3)/lz;
xy1(:,4) = xy1(:,4)/lx4;



%%
X1            =           xy1(:,1);
NxNy1         =           length(X1);
X2            =           xy_s(:,1);
NxNy2         =           length(X2);
weight        =           zeros(NxNy1,k);

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

    w = RBF_interp_PHS_pol_weights_4d (xk1, xe1, m,d);      % RBF-PHS method

    weight(m1,1:k0) = w;

    %
end


Idx_x                                  =       reshape((1:NxNy1).*ones(k,1),1,[]);
Idx_y                                  =       reshape(Nearest_Idx_interp(:,1:k)',1,[]);
    W                                  =       reshape(squeeze(weight(:,:))',1,[]);
    D0                                 =       sparse(Idx_x,Idx_y,W,NxNy1,NxNy2);


end