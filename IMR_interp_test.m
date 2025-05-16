
addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')
addpath('functions/RBF')


%%  2D parameter space (NHKV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%Sampling from Prior%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_true = 'neoHook';  

mu_theta          =  [15090 0.209 0];
sigma_theta       =  ([4350 10 0; 10 0.1 0; 0 0 0]).^2;

N_target = 5000;

theta_target =  mvnrnd(mu_theta,sigma_theta,N_target);


%%%%%%%%%%%%%%%%%%%%%%%%%%%Sampling for Batch Samples%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


boundary_idx = (theta_target(:,1) >mu_theta(1)+sqrt(sigma_theta(1,1))*2.5| theta_target(:,1) <mu_theta(1)-sqrt(sigma_theta(1,1))*2.5) | (theta_target(:,2) >mu_theta(2)+sqrt(sigma_theta(2,2))*2.5| theta_target(:,2) <mu_theta(2)-sqrt(sigma_theta(2,2))*2.5);

idx_all = 1:N_target;

params_idx = [randi([1 N_target],400,1); (idx_all(boundary_idx ))'];

params_idx  = unique( params_idx,'rows');

theta_params =  theta_target(params_idx,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%Simulations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

Design = [480, 0.15];

  [t,y2,U2,dt] = IMR_simulations(theta_target,Design,model_true,60);

toc

figure; 
plot(t,y2,'k-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%Seach Nearest Neighibors%%%%%%%%%%%%%%%%%%%%%%%%%

k= 40;

[Nearest_Idx] = nearest_interp(theta_target(:,1:2),theta_params(:,1:2),k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%RBF Interpolation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I]     =     RBF_PHS_interp(theta_target(:,1:2),theta_params(:,1:2),Nearest_Idx,20,'1',3,2);


y_interp = I*y2(params_idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%Test for target parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_idx = 8723;


[t0,y_target,~,~] = IMR_simulations(theta_target(target_idx,:),Design,model_true,60);


sigma_w = @(y) noise_fun(t0,y);

figure; 

subplot(2,1,1)

plot(t0,y_target,'ro-'); hold on;

plot(t,y_interp(target_idx,:),'b*-'); hold on;

% plot(t,y_interp_2(target_idx,:),'m.-'); hold on;

% plot(t,y_interp1(target_idx,:),'m.-'); hold on;

plot(t0,y_target+sqrt(sigma_w(y_target))*2,'r--');
plot(t0,y_target-sqrt(sigma_w(y_target))*2,'r--');

 plot(t,y2(Nearest_Idx(target_idx,1:20),:),'k:'); hold on;

ylim([0 1])

legend('target','interpolation')

xlabel('t')

ylabel('R')


subplot(2,1,2)

scatter(theta_params(:,1),theta_params(:,2),'b.'); axis square tight; hold on; 

scatter(theta_target(target_idx,1), theta_target(target_idx,2),'r*');

scatter(theta_params(Nearest_Idx(target_idx,1:20),1),theta_params(Nearest_Idx(target_idx,1:20),2),'ko');

xlabel('G');
ylabel('\mu')

legend('samples','target','neighbors')


       Sigma = diag(sigma_w(y_target));

     (norm_fun_log(y_target, Sigma, y_interp(target_idx,:)))


       Sigma = diag(sigma_w(mean(y2,1)));

     (norm_fun_log(mean(y2,1), Sigma, mean(y_interp,1)))


scatter(theta_params(:,1),theta_params(:,2),'r*'); axis square tight; hold on; 

scatter(theta_target(:,1), theta_target(:,2),'b.');

xlabel('G');
ylabel('\mu')

legend('samples','target','neighbors')


%%  3D parameter space (qKV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%Sampling from Prior%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_true    = 'fung';  
mu_theta      =  [2770 0.286 0.28];
sigma_theta   =  ([300 0 0; 0 0.106 0; 0 0 0.1]).^2;
N_target      =  500;
theta_target  =  mvnrnd(mu_theta,sigma_theta,N_target);


%%%%%%%%%%%%%%%%%%%%%%%%%%%Sampling for Batch Samples%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundary_idx = (theta_target(:,1) >mu_theta(1)+sqrt(sigma_theta(1,1))*2| theta_target(:,1) <mu_theta(1)-sqrt(sigma_theta(1,1))*2)...
     | (theta_target(:,2) >mu_theta(2)+sqrt(sigma_theta(2,2))*2| theta_target(:,2) <mu_theta(2)-sqrt(sigma_theta(2,2))*2)...
      | (theta_target(:,3) >mu_theta(3)+sqrt(sigma_theta(3,3))*2| theta_target(:,3) <mu_theta(3)-sqrt(sigma_theta(3,3))*2) ;


idx_all = 1:N_target;

params_idx = [randi([1 N_target],200,1); (idx_all(boundary_idx ))'];

params_idx  = unique( params_idx,'rows');

 theta_params =  theta_target(params_idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%Simulations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

Design = [480, 0.15];

  [t,y2,U2,dt] = IMR_simulations(theta_target,Design,model_true,60);

  toc

figure; 

plot(t,y2,'k-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%Seach Nearest Neighibors%%%%%%%%%%%%%%%%%%%%%%%%%

k= 40;

[Nearest_Idx] = nearest_interp_3d(theta_target(:,:),theta_params(:,:),k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%RBF Interpolation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I]     =     RBF_PHS_interp_3d(theta_target(:,:),theta_params(:,:),Nearest_Idx,15,'1',3,2);

y_interp = I*y2(params_idx,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%Test for target parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_idx = 5017;

[t,y_target,U_target,dt] = IMR_simulations(theta_target(target_idx,:),Design,model_true,60);


figure;

subplot(2,1,1)

plot(t,y_target,'ro-'); hold on;

plot(t,y_interp(target_idx,:),'b*-'); hold on;

plot(t,y2(Nearest_Idx(target_idx,1:15),:),'k:'); hold on;

ylim([0 1])

legend('target','interpolation')

xlabel('t')

ylabel('R')


subplot(2,1,2)

scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,3),'b.'); axis square tight; hold on;

scatter3(theta_target(target_idx,1), theta_target(target_idx,2),theta_target(target_idx,3),'r*');

% scatter(theta_target(:,1), theta_target(:,2),'r*');

scatter3(theta_params(Nearest_Idx(target_idx,1:15),1),theta_params(Nearest_Idx(target_idx,1:15),2),theta_params(Nearest_Idx(target_idx,1:15),3),'ko');

xlabel('G');
ylabel('\mu')
zlabel('\alpha')

legend('samples','target','neighbors')



%%  4D

model_true = 'neoHook';  

mu_theta          =  [15090 0.209 0];
sigma_theta       =  ([4350 0 0; 0 0.1 0; 0 0 0]).^2;

N_target = 2000;

theta_target =  mvnrnd(mu_theta,sigma_theta,N_target);
design_target = [unifrnd(100,1000,[N_target 1]) unifrnd(0.14,0.3,[N_target 1])];


boundary_idx = (theta_target(:,1) >mu_theta(1)+sqrt(sigma_theta(1,1))*2| theta_target(:,1) <mu_theta(1)-sqrt(sigma_theta(1,1))*2) | (theta_target(:,2) >mu_theta(2)+sqrt(sigma_theta(2,2))*2| theta_target(:,2) <mu_theta(2)-sqrt(sigma_theta(2,2))*2);

idx_all = 1:N_target;

params_idx = [(1:900)'; (idx_all(boundary_idx ))'];

params_idx  = unique( params_idx,'rows');

 theta_params =  theta_target(params_idx,:);
 design_params =  design_target(params_idx,:);

Y = zeros(N_target,61);

tic

parfor j = 1:N_target

   disp(['simulation #:' num2str(j)])

  [t,y2,U2,dt] = IMR_simulations(theta_target(j,:),design_target(j,:),model_true,60);

  Y(j,:) = y2;

end

toc

t = linspace(0,3,61);


figure; 

plot(t,Y,'k-')


k= 40;

[Nearest_Idx] = nearest_interp([theta_target(:,1:2), design_target],[theta_params(:,1:2), design_params],k);


[I]     =     RBF_PHS_interp_4d([theta_target(:,1:2), design_target],[theta_params(:,1:2), design_params],Nearest_Idx,40,'1',3,2);


y_interp = I*Y(params_idx,:);


target_idx = 1932;

[t,y_target,U_target,dt] = IMR_simulations(theta_target(target_idx,:),design_target(target_idx,:),model_true,60);


figure; 

subplot(2,1,1)


plot(t,y_target,'ro-'); hold on;

plot(t,y_interp(target_idx,:),'b*-'); hold on;

 plot(t,Y(Nearest_Idx(target_idx,1:30),:),'k:'); hold on;

ylim([0 1])

legend('target','interpolation')

xlabel('t')

ylabel('R')



subplot(2,1,2)

scatter3(theta_params(:,1),theta_params(:,2),design_params(:,1),10, design_params(:,2),'+'); axis square tight; hold on; 
hcb=colorbar;
hcb.Title.String = "Req";

caxis([0.14 0.3] )


scatter3(theta_target(target_idx,1), theta_target(target_idx,2),design_target(target_idx,1),50,  design_target(target_idx,2),'*'); axis square tight; hold on; 


 % scatter(theta_target(:,1), theta_target(:,2),'r*');

scatter3(theta_params(Nearest_Idx(target_idx,1:30),1),theta_params(Nearest_Idx(target_idx,1:30),2),design_params(Nearest_Idx(target_idx,1:30),1),50, design_params(Nearest_Idx(target_idx,1:30),2),'o'); axis square tight; hold on; 

xlabel('G');
ylabel('\mu')
zlabel('We')


legend('samples','target','neighbors')






%%  error analysis


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


model_true = 'neoHook';  

mu_theta          =  [15090 0.209 0];
sigma_theta       =  ([4350 0 0; 0 0.1 0; 0 0 0]).^2;

% design_available = [894 0.21;876 0.25; 875 0.21; 833 0.24; 814 0.27; 789 0.16; 776 0.17; 757 0.15; 619 0.30; 600 0.30; 539 0.16;...
%     526 0.18; 516 0.22; 477 0.16; 446 0.30; 395 0.25; 339 0.19; 188 0.25; 187 0.24; 178 0.18; 150 0.30];


design_available = [203 0.27];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(2*Sigma))*(x-mu)' );  % removed the constant here



figure;

cmap = colormap(jet);

for j = 1:length(design_available)

    % j = 8;

    Design = design_available(j,:);

    data_filename =  ['data_10k/Rt', '_We',num2str(round(Design(1)),'%i'),'_Req', ...
        num2str(Design(2),'%.2f') '.mat']; % name of file containing R vs T data

    load(data_filename);

    y2 = data{1}; theta_target = data{2};

    N_target = length(theta_target);

     k_std = 3;

    % boundary_idx = (theta_target(:,1) >mu_theta(1)+sqrt(sigma_theta(1,1))*k_std| theta_target(:,1) <mu_theta(1)-sqrt(sigma_theta(1,1))*k_std) | (theta_target(:,2) >mu_theta(2)+sqrt(sigma_theta(2,2))*k_std| theta_target(:,2) <mu_theta(2)-sqrt(sigma_theta(2,2))*k_std);

     boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;


    idx_all = 1:N_target;

    error = zeros(N_target,20);
    error2 = zeros(N_target,20);

    % 
    % [~,y_center,U_center,~] = IMR_simulations(mu_theta,Design,model_true,60);


    for i = 1:150


        params_idx = [randi([1 N_target],ceil(((i)*1)/4)^2+10,1); (idx_all(boundary_idx ))'];

        params_idx  = unique( params_idx,'rows');

        theta_params =  theta_target(params_idx,:);

        theta_params_rnd =  [unifrnd(min(theta_target(:,1)),max(theta_target(:,1)),length(theta_params),1)...
                         unifrnd(min(theta_target(:,2)),max(theta_target(:,2)),length(theta_params),1) ...
                         unifrnd(min(theta_target(:,3)),max(theta_target(:,3)),length(theta_params),1)]; 

        k= 20;

        [Nearest_Idx] = nearest_interp(theta_target(:,1:2),theta_params(:,1:2),k);

        [I]     =     RBF_PHS_interp(theta_target(:,1:2),theta_params(:,1:2),Nearest_Idx,20,'1',3,2);

        y_interp = I*y2(params_idx,:);

        y_interp(params_idx,:) = y2(params_idx,:);


        Sigma_mean = 0;
        y_mean = 0;

        for k = 1:N_target
            Sigma = diag(sigma_w(y2(k,:)));
            Sigma2 = diag(sigma_w(y_interp(k,:)));
            Sigma_mean = Sigma_mean +(Sigma+(y2(k,:)-mean(y2,1))'*(y2(k,:)-mean(y2,1)))/N_target;
            % Sigma = eye(61,61);
            error(k,i)   = (norm_fun_log(y2(k,:), Sigma, y_interp(k,:)));
            error_2(k,i) = (norm_fun_log(y2(k,:), Sigma2, y_interp(k,:)))+...
                 1/2*(trace((Sigma2)\Sigma)-length(y2(k,:))+log(det(Sigma2)/det(Sigma)));

        end


        for k = 1:100:N_target
             error_2(k,i) = error(k,i)-(norm_fun_log(mean(y2,1), Sigma_mean,mean(y_interp)))/1^2;
            % error_2(k,i) = error_2(k,i)-log(mean(exp(error_2(:,i)),1));

            scatter(k/N_target, mean((error(1:k,i)) ),[],(j-1),'.'); hold on; drawnow;
        end


       %  scatter(length(theta_params)/N_target, sqrt(sum(error(:,i))/N_target),'r.'); hold on;

         scatter(length(theta_params)/N_target, mean((error(:,i)) ),[],(j-1),'.'); hold on;

         if mean((error_2(:,i)) )>0
          scatter(length(theta_params)/N_target, mean((error_2(:,i)) ),[],(j-1),'r+'); hold on;
         end

        axis square tight;



        % figure;
        % 
        % for k = 1:N_target
        % 
        %     plot(k, mean((error(1:k,i)) ),'k.'); hold on; 
        % end

        % xlim([1e-2 1]);

        % ylim([0.14 0.3]);
        %
        % zlim([0 3]);

        xscale('log');
        yscale('log');

        drawnow;

    end

end

xlabel('$N_B/N_S$','Interpreter','latex')
ylabel('$\sqrt{(\delta Y)^T\Sigma^{-1}(\delta Y)}$','Interpreter','latex')




%%  2D design

model_true = 'neoHook';  

mu_theta          =  [15090 0.209 0];

N_target = 2000;

design_target = [unifrnd(100,1000,[N_target 1]) unifrnd(0.14,0.3,[N_target 1])];

% theta_target =  mvnrnd(mu_theta,sigma_theta,N_target);
% 
% boundary_idx = (theta_target(:,1) >mu_theta(1)+sqrt(sigma_theta(1,1))*2| theta_target(:,1) <mu_theta(1)-sqrt(sigma_theta(1,1))*2) | (theta_target(:,2) >mu_theta(2)+sqrt(sigma_theta(2,2))*2| theta_target(:,2) <mu_theta(2)-sqrt(sigma_theta(2,2))*2);
% 
% idx_all = 1:N_target;
% 
% params_idx = [(1:200)'; (idx_all(boundary_idx ))'];
% 
% params_idx  = unique( params_idx,'rows');


params_idx = 1:200;

 design_params =  design_target(params_idx,:);

Y = zeros(N_target,61);

tic

parfor j = 1:N_target

  [t,y2,U2,dt] = IMR_simulations(mu_theta,design_target(j,:),model_true,60);

  Y(j,:) = y2;

end

  toc

 t  = linspace(0,3,61);

figure; 

plot(t,Y,'k-')


k= 40;

[Nearest_Idx] = nearest_interp(design_target(:,1:2),design_params(:,1:2),k);


[I]     =     RBF_PHS_interp(design_target(:,1:2),design_params(:,1:2),Nearest_Idx,15,'1',3,2);

[Dx]     =     RBF_PHS_interp(design_target(:,1:2),design_params(:,1:2),Nearest_Idx,20,'x',3,3);

[Dy]     =     RBF_PHS_interp(design_target(:,1:2),design_params(:,1:2),Nearest_Idx,20,'y',3,3);


y_interp = I*Y(params_idx,:);

dydd1    = Dx*Y(params_idx,:);
dydd2    = Dy*Y(params_idx,:);

target_idx = 913;

[t,y_target,U_target,dt] = IMR_simulations(mu_theta,design_target(target_idx,:),model_true,60);


figure; 

subplot(2,1,1)


plot(t,y_target,'ro-'); hold on;

plot(t,y_interp(target_idx,:),'b*-'); hold on;

 plot(t,Y(Nearest_Idx(target_idx,1:15),:),'k:'); hold on;

ylim([0 1])

legend('target','interpolation')

xlabel('t')

ylabel('R')



subplot(2,1,2)

scatter(design_params(:,1),design_params(:,2),'b.'); axis square tight; hold on; 

scatter(design_target(target_idx,1), design_target(target_idx,2),'r*');

 % scatter(theta_target(:,1), theta_target(:,2),'r*');

scatter(design_params(Nearest_Idx(target_idx,1:15),1),design_params(Nearest_Idx(target_idx,1:15),2),'ko');

xlabel('We');
ylabel('R_{\infty}')

legend('samples','target','neighbors')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dd1 =1e-2; dd2 = 1e-5;


target_idx = 913;

[~,y_target_P,~,~] = IMR_simulations(mu_theta,[design_target(target_idx,1)+dd1, design_target(target_idx,2)],model_true,60);

[~,y_target_M,~,~] = IMR_simulations(mu_theta,[design_target(target_idx,1)-dd1, design_target(target_idx,2)],model_true,60);


dydd1_target = (y_target_P-y_target_M)/2/dd1;

[~,y_target_P,~,~] = IMR_simulations(mu_theta,[design_target(target_idx,1), design_target(target_idx,2)+dd2],model_true,60);

[~,y_target_M,~,~] = IMR_simulations(mu_theta,[design_target(target_idx,1), design_target(target_idx,2)-dd2],model_true,60);


dydd2_target = (y_target_P-y_target_M)/2/dd2;



figure; 

subplot(3,1,1)


plot(t,dydd1_target,'ro-'); hold on;

plot(t,dydd1(target_idx,:),'b*-'); hold on;

 % plot(t,Y(Nearest_Idx(target_idx,1:15),:),'k:'); hold on;

% ylim([0 1])

legend('target','RBF-FD')

ylabel('\partial_{d1} R')


subplot(3,1,2)


plot(t,dydd2_target,'ro-'); hold on;

plot(t,dydd2(target_idx,:),'b*-'); hold on;

 % plot(t,Y(Nearest_Idx(target_idx,1:15),:),'k:'); hold on;

% ylim([0 1])

legend('target','RBF-FD')

xlabel('t')

ylabel('\partial_{d2} R')




subplot(3,1,3)

scatter(design_params(:,1),design_params(:,2),'b.'); axis square tight; hold on; 

scatter(design_target(target_idx,1), design_target(target_idx,2),'r*');

 % scatter(theta_target(:,1), theta_target(:,2),'r*');

scatter(design_params(Nearest_Idx(target_idx,1:15),1),design_params(Nearest_Idx(target_idx,1:15),2),'ko');

xlabel('We');
ylabel('R_{\infty}')

legend('samples','target','neighbors')





