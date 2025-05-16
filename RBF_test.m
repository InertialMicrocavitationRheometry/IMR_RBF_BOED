
k= 40;

[Nearest_Idx] = nearest_interp(xy1,xy1,k,10);


% 
    [D]   =   RBF_PHS_interp_all(xy1,xy1,Nearest_Idx, 36, 3, 3);

%    
% 
% [D]   =   RBF_PHS_interp_all(xy1,xy1,Nearest_Idx, 17, 3, 2);

Dx    =   D{1};
Dy    =   D{2};


L0    =   D{3};


Dxx   =   D{4};
Dyy   =   D{5};



 Nearest_Idx_nb =  xy1(:,1)< x_min+.2 | xy1(:,2)< y_min+0.8 | xy1(:,2)> y_max-0.8 ;

Idx = (1:length(xy1))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);


xy_nb = xy1(Nearest_Idx_nb,:);

[Dx_nb]     =     RBF_PHS_interp(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),25,'x',3,2);
[Dy_nb]     =     RBF_PHS_interp(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),25,'y',3,2);
[Dxx_nb]    =     RBF_PHS_interp(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),25,'xx',3,2);
[Dyy_nb]    =     RBF_PHS_interp(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),25,'yy',3,2);


Dx(Nearest_Idx_nb,:)   = Dx_nb;
Dy(Nearest_Idx_nb,:)   = Dy_nb;
Dxx(Nearest_Idx_nb,:)  = Dxx_nb;
Dyy(Nearest_Idx_nb,:)  = Dyy_nb;


