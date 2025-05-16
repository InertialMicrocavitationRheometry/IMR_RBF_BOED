function [Nearest_Idx] = nearest_interp_3d(xy,xy_s,k)

% for each node in xy, find its k nearest nerighbor in xy_s


X1            =           xy(:,1);

NxNy1         =           length(X1);




%% standardize the searching region

lx = max(xy(:,1))- min(xy(:,1));

ly = max(xy(:,2))- min(xy(:,2));

lz = max(xy(:,3))- min(xy(:,3));


if lx<1e-6
    lx=1e-6;
end

 if ly<1e-6
    ly=1e-6;   
 end

  if lz<1e-6
    lz=1e-6;   
 end


 lx = 1;
 ly = 1;
lz = 1;

     % lz = 10;

xy_s(:,1) = xy_s(:,1)/lx;
xy_s(:,2) = xy_s(:,2)/ly;
xy_s(:,3) = xy_s(:,3)/lz;

xy(:,1) = xy(:,1)/lx;
xy(:,2) = xy(:,2)/ly;
xy(:,3) = xy(:,3)/lz;




%%  Find nearest nodes

Nearest_Idx     =       zeros(NxNy1,k);



disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Find nearest ' num2str(k) ' nodes'])


for j=1:NxNy1

    if mod(j,100)==0
        disp(['j= ' num2str(j)])
    end

    Node_j                       =       xy(j,:);


    r2_j                         =       (xy_s-Node_j).^2;
    r2_j                         =       r2_j(:,1)+r2_j(:,2)+r2_j(:,3);



    % Neighbor region

    Idx_n                        =       find(r2_j<0.3^2);

    if length(Idx_n)<k


        Idx_n                        =       find(r2_j<30000^2);

    end




    % search k nearest nodes in the neighbor region

    Idx_0                        =       knnsearch(Node_j,xy_s(Idx_n,:),k);


    Idx                          =       Idx_n(Idx_0(1:k));



    Nearest_Idx(j,:)             =       Idx;



end


end
