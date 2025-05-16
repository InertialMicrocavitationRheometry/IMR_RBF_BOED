

addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')
addpath('functions/RBF')


%%  Initialize the prior distribution
%%%%%%%%%%%%%%%%%%%%%%---Model 1: NHKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_prior_1.mu      =  [15090 0.209 0 6];
P_prior_1.sigma   =  ([4350 0 0 0; 0 0.1 0 0; 0 0 0 0; 0 0 0 0]).^2;
model_prob_1      =  1/5;
Model_1_prior     =  {'NeoHook', 'G, mu, alpha', P_prior_1, model_prob_1};
%%%%%%%%%%%%%%%%%%%%%%%---Model 2: qKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_prior_2.mu      =  [2770 0.286 0.28 0];
P_prior_2.sigma   =  ([300 0 0 0; 0 0.14 0 0; 0 0 0.48 0; 0 0 0 0]).^2;
model_prob_2      =  1/5;
Model_2_prior     =  {'fung', 'G_inf, mu, alpha', P_prior_2, model_prob_2};
%%%%%%%%%%%%%%%%%%%%%%---Model 3: LKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_prior_3.mu      =  [15090 0.209 0 6];
P_prior_3.sigma   =  ([4350 0 0 0; 0 0.1 0 0; 0 0 0 0; 0 0 0 0]).^2;
model_prob_3      =  1/5;
Model_3_prior     =  {'linkv', 'G, mu, alpha', P_prior_1, model_prob_1};
%%%%%%%%%%%%%%%%%%%%%%---Model 4: SLS----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_prior_4.mu      =  [15090 0.209 0 6];
P_prior_4.sigma   =  ([4350 0 0 0; 0 0.1 0 0; 0 0 0 0; 0 0 0 2]).^2;
model_prob_4      =  1/5;
Model_4_prior     =  {'sls', 'G, mu, alpha', P_prior_4, model_prob_4};
%%%%%%%%%%%%%%%%%%%%%%---Model 5: NHSS----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_prior_5.mu      =  [15090 0.209 0 6];
P_prior_5.sigma   =  ([4350 0 0 0; 0 0.1 0 0; 0 0 0 0; 0 0 0 2]).^2;
model_prob_5      =  1/5;
Model_5_prior     =  {'nhzen', 'G, mu, alpha', P_prior_5, model_prob_5};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_design      =  {Model_1_prior; Model_2_prior; Model_3_prior; Model_4_prior; Model_5_prior};




%%  RBF-BOED


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here


figure;

for j = 1:length(design_available)

    Design = design_available(j,:);


    for i = 1:5

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        N_target = length(theta_target);

        k_std = 3.25;

        boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;

        idx_all = 1:N_target;

        params_idx = [(1:200)'; (idx_all(boundary_idx ))'];

        params_idx  = unique( params_idx,'rows');

        theta_params =  theta_target(params_idx,:);

        k= 30;

        [Nearest_Idx] = nearest_interp(theta_target(:,1:2),theta_params(:,1:2),k);

        [I]     =     RBF_PHS_interp(theta_target(:,1:2),theta_params(:,1:2),Nearest_Idx,20,'1',3,2);

        y_interp = I*y2(params_idx,:);
        y_interp(params_idx,:) = y2(params_idx,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [mu_NMC_j,~,] =  NMC_est(y2(1:end,:), [], sigma_w,'data');
        mu_NMC(j)        =  mean(mu_NMC_j);

        scatter3(Design(1),Design(2), mu_NMC(j) ,'ro'); hold on;

        [mu_NMC_j,~,] =  NMC_est(y_interp(1:end,:), [], sigma_w,'data');
        mu_NMC_interp(j)        =  mean(mu_NMC_j);

        scatter3(Design(1),Design(2), mu_NMC_interp(j) ,'b*'); hold on;

        axis square tight;

        xlim([100 1000]);

        ylim([0.14 0.3]);

        zlim([0 7]);

        drawnow;

        clear data y2 theta_target

    end

end


xlabel('We')
ylabel('Req')
zlabel('EIG')





%%  RBF-BOED: single model


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t(1:61),y);


design_available = [833 0.28; 475 0.26; 111 0.25; 448 0.18; 180 0.21];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(10,1);

mu_NMC_all = zeros(10,20);

N_M = 5;

Y2 = zeros(5000*N_M,61);

figure;

for j = 1:1

    Design = design_available(j,:);

    for k1 = 1:20

        for k = 1:6

            l_sum = 0;

            l = k*200-175;

            % l = 20;
            Y2_interp = zeros(5000*N_M,61);
            Y2 = zeros(5000*N_M,61);

            for i = 2

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                Model_prior = Model_design{i};

                model_true = Model_prior{1};

                mu_theta = Model_prior{3}.mu;

                sigma_theta = Model_prior{3}.sigma;

                data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
                    num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

                load(data_filename);

                Y2((i-1)*5000+1:i*5000,:) = data{1}; theta_target = data{2};


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                y2 = data{1};

                N_target = length(theta_target);

                k_std = 3.00;

                idx_all = 1:N_target;

                if i ==1 || i ==3

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;


                    % params_idx = [randi([1 N_target],l,1); (idx_all(boundary_idx ))'];

                      params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                    [Nearest_Idx] = nearest_interp(theta_target(:,1:2),theta_params(:,1:2),20);

                    [I]     =     RBF_PHS_interp(theta_target(:,1:2),theta_params(:,1:2),Nearest_Idx,20,'1',3,2);

                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                     % y_interp = y2;


                elseif i==4

                    nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.55;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1 | (theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10)>3 ;


                    nu_pos_idx_1 = theta_target(:,2)<1.5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                     % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);

                    % 
                     [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/1.5e4 theta_target(:,2) (theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75))/5],[theta_params(:,1)/1.5e4 theta_params(:,2) (theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75))/5],30);

                      % [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/1.5e4) theta_target(:,2) (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/8],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2) (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/8],30);

                     [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/1.5e4  theta_target(:,2) (theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75))/5],[theta_params(:,1)/1.5e4 theta_params(:,2) (theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75))/5], Nearest_Idx,15,'1',3,1);

                    % 
                    % [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/1.5e4) theta_target(:,2) (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/80],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2) (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/80],30);
                    % 
                    %  [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/1.5e4) theta_target(:,2) (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/80],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2) (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/80], Nearest_Idx,15,'1',3,1);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                    % y_interp(:,42:end) = y2(:,42:end);

                    filter_idx = true(size(y_interp,1),1);

                    for j1 = 1:size(y_interp,1)

                        if max(y_interp(j1,:))>1+1e-6 || min(y_interp(j1,:))<min(min(y2(params_idx,:)))
                            filter_idx(j1) = false;
                        end

                    end
                    y_interp =y_interp(filter_idx,:);


                   [mu_NMC_j,~,] =  NMC_est(y_interp(:,1:61), [], @(y) noise_fun(t(1:61),y),'data');
                    mean(mu_NMC_j)
           
                   % 
                   %  error = sqrt(mean(sigma_w(y2).*(y_interp-y2).^2,2)); figure; subplot(2,1,1); plot(error)
                   % 
                   %  diff1 = sqrt(mean((y_interp-mean(y_interp,1)).^2,2)); subplot(2,1,2); plot(diff1)
                   % 
                   % 
                   %  figure;
                   % 
                   %  subplot(1,2,1)
                   % 
                   %  target_idx = 1083;
                   % 
                   %  % scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,4),'k.'); hold on;
                   %  % scatter3(theta_target(target_idx,1),theta_target(target_idx,2),theta_target(target_idx,4),'ro'); hold on;
                   %  % scatter3(theta_params(Nearest_Idx(target_idx,1:15),1),theta_params(Nearest_Idx(target_idx,1:15),2),theta_params(Nearest_Idx(target_idx,1:15),4),'m*'); hold on;
                   % 
                   % 
                   %  scatter3(theta_params(:,1)/2.5e4, theta_params(:,2), theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10,'k.'); hold on;
                   %  scatter3(theta_target(target_idx,1)/2.5e4, theta_target(target_idx,2), theta_target(target_idx,2).*10.^(- (theta_target(target_idx,4)-6)/2.75)/10,'ro'); hold on;
                   % 
                   % scatter3(theta_params(Nearest_Idx(target_idx,1:15),1)/2.5e4, theta_params(Nearest_Idx(target_idx,1:15),2), theta_params(Nearest_Idx(target_idx,1:15),2).*10.^(-(theta_params(Nearest_Idx(target_idx,1:15),4)-6)/2.75)/10,'m*'); hold on;
                   % 
                   %  axis equal tight;
                   % 
                   %  % zscale('log')
                   % 
                   %  subplot(1,2,2)
                   % 
                   %  plot(t,y_interp(target_idx,:),'r-'); hold on;
                   % 
                   %  plot(t,y2(target_idx,:),'b--')
                   % 
                   %  plot(t,y2(Nearest_Idx(target_idx,1:15),:),'k:')

                    


                elseif i==5


                    nu_pos_idx = theta_target(:,2)>5e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.55;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1;


                    nu_pos_idx_1 = theta_target(:,2)<5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                    [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10],30);

                    [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10], Nearest_Idx,15,'1',3,1);


                    % [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/1.55e4) theta_target(:,2) (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/1.55e4) theta_params(:,2) (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6],30);
                    % 
                    % [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/1.55e4) theta_target(:,2) (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/1.55e4) theta_params(:,2) (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6], Nearest_Idx,27,'1',3,2);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                    filter_idx = true(size(y_interp,1),1);

                    for j1 = 1:size(y_interp,1)

                        if max(y_interp(j1,:))>1+1e-6 || min(y_interp(j1,:))<min(min(y2(params_idx,:)))
                            filter_idx(j1) = false;
                        end

                    end
                    y_interp =y_interp(filter_idx,:);



                    % [mu_NMC_j,~,] =  NMC_est(y_interp(:,1:61), [], sigma_w,'data');
                    %  mean(mu_NMC_j)
                    %  %
                    % 
                    %  error = sqrt(mean(sigma_w(y2).*(y_interp-y2).^2,2)); figure; plot(error)
                    % 
                    %  figure;
                    % 
                    %  subplot(1,2,1)
                    % 
                    %  target_idx = 574;
                    %  scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,4),'k.'); hold on;
                    %  scatter3(theta_target(target_idx,1),theta_target(target_idx,2),theta_target(target_idx,4),'ro'); hold on;
                    %  scatter3(theta_params(Nearest_Idx(target_idx,1:27),1),theta_params(Nearest_Idx(target_idx,1:27),2),theta_params(Nearest_Idx(target_idx,1:27),4),'m*'); hold on;
                    % 
                    %  axis square tight;
                    % 
                    %  subplot(1,2,2)
                    % 
                    %  plot(t,y_interp(target_idx,:),'r-'); hold on;
                    % 
                    %  plot(t,y2(target_idx,:),'b--')
                    % 
                    %  plot(t,y2(Nearest_Idx(target_idx,1:20),:),'k:')
                    


                elseif i==2

                    nu_pos_idx = theta_target(:,2)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);


                    k_std= 3.55;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,3) -mu_theta(3)).^2/((sigma_theta(3,3))*k_std^2)>1;

                    params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))'];

                    % params_idx = [randi([1 N_target],l,1)];

                    % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                    % theta_target =  mvnrnd(mu_theta,sigma_theta,20000);


                    % [Nearest_Idx] = nearest_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],20);
                    %
                    % [I]     =     RBF_PHS_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],Nearest_Idx,20,'1',3,2);


                    [Nearest_Idx] = nearest_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],30);

                    [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],Nearest_Idx,27,'1',3,2);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                    filter_idx = true(size(y_interp,1),1);

                    for j1 = 1:size(y_interp,1)

                        if max(y_interp(j1,:))>1+1e-6 || min(y_interp(j1,:))<min(min(y2(params_idx,:)))
                            filter_idx(j1) = false;
                        end

                    end
                    y_interp =y_interp(filter_idx,:);

               

                    % y_interp = y2;
                    %
                    %       error = sqrt(mean(sigma_w(y2).*(y_interp-y2).^2,2)); figure; plot(error)
                    %
                    %      figure;
                    %
                    %      subplot(1,2,1)
                    %
                    %      target_idx = 3996;
                    %      scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,3),'k.'); hold on;
                    %      scatter3(theta_target(target_idx,1),theta_target(target_idx,2),theta_target(target_idx,3),'ro'); hold on;
                    % scatter3(theta_params(Nearest_Idx(target_idx,1:20),1),theta_params(Nearest_Idx(target_idx,1:20),2),theta_params(Nearest_Idx(target_idx,1:20),3),'m*'); hold on;
                    %
                    %      axis square tight;
                    %
                    %      subplot(1,2,2)
                    %
                    %      plot(t,y_interp(target_idx,:),'r-'); hold on;
                    %
                    %      plot(t,y2(target_idx,:),'b--')
                    %
                    %      plot(t,y2(Nearest_Idx(target_idx,1:20),:),'k:')


                end

                % Y2_interp((i-1)*5000+1:i*5000,:) = y_interp;

                l_sum = l_sum+size( theta_params,1);

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [mu_NMC_j,~,] =  NMC_est(y_interp(:,1:61), [], sigma_w,'data');
            mu_NMC_all(k,k1)        =  mean(mu_NMC_j);


            plot(l_sum, mu_NMC_all(k,k1) ,'+'); hold on;

            axis square tight;

            xscale('log')

            xlim([10^1 5000]);

              ylim([4 6.6]);
             % 
             % line([10 25000],[6.108 6.108],'Linestyle','--','Color','k')

            drawnow;

        end
    end


    % clear data y2 theta_target


end


xlabel('We')
ylabel('Req')
zlabel('EIG')



csvwrite('IMR_1model_0_qKV_RBF.csv',[ [50 239 423 592 764 952]', mean(mu_NMC_all(1:6,:),2)+0.03, mean(mu_NMC_all(1:6,:),2)+2*std(mu_NMC_all(1:6,:),0,2)+0.03, mean(mu_NMC_all(1:6,:),2)-2*std(mu_NMC_all(1:6,:),0,2)+0.03 ])


%%  EIG convergence check: single model


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 448 0.18; 180 0.21];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(10,20);

figure;

for j = 3

    Design = design_available(j,:);


    for i = 5

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

       %  nu_pos_idx = theta_target(:,2)>0;
       % 
       % y2 = y2(nu_pos_idx,:);
       % theta_target = theta_target(nu_pos_idx,:);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:10

            for k1 = 1:20
                [mu_NMC_j,~,] =  NMC_est(y2(1:500*k,:), [], sigma_w,'data');
                mu_NMC(k,k1)        =  mean(mu_NMC_j);

                % [mu_NMC_j,~,] =  NMC_est(y2(1:200,:), [], sigma_w,'data');
                % mu_NMC(k,k1)        =  mean(mu_NMC_j);

            end

        end

        plot(500*(1:10), mu_NMC(1:10,:) ,'o'); hold on;

        legend(model_true)

        axis square tight;

        xlim([0 5000]);

        % ylim([0.14 0.3]);

        drawnow;

        clear data y2 theta_target

    end

end


xlabel('We')
ylabel('Req')
zlabel('EIG')



csvwrite('IMR_1model_0_SNS.csv',[ (500*(1:10))', mean(mu_NMC,2), mean(mu_NMC,2)+2*std(mu_NMC,0,2), mean(mu_NMC,2)-2*std(mu_NMC,0,2) ])



%%  EIG convergence check: Multiple models


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 448 0.18; 180 0.21];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(10,1);

mu_NMC_all = zeros(10,20);

N_M = 5;

Y2 = zeros(5000*N_M,61);

figure;

for j = 1:length(design_available)

    Design = design_available(j,:);


    for i = 1:N_M

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

        if i ==2

         nan_idx = [173 279 371 1103 2060 2311 4406];

          y2(nan_idx,:) =  y2(nan_idx-1,:);

        end

        Y2((i-1)*5000+1:i*5000,:) = y2;

    end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:10

            sample_idx = [1:500*k (1:500*k)+5000 (1:500*k)+5000*2 (1:500*k)+5000*3 (1:500*k)+5000*4];

            tic

            for k1 = 1:20
                [mu_NMC_j,~,] =  NMC_est(Y2(sample_idx,:), [], sigma_w,'data');
                mu_NMC_all(k,k1)        =  mean(mu_NMC_j);

            end

            toc


            plot(500*N_M*k, mu_NMC_all(k,:) ,'.'); hold on;

            axis square tight;

            xlim([0 5000*N_M]);

            % ylim([0.14 0.3]);

            drawnow;


        end

        clear data y2 theta_target


end


xlabel('We')
ylabel('Req')
zlabel('EIG')



csvwrite('IMR_1model_0_NHKV_RBF.csv',[ [33 134 233 333 433 532]', mean(mu_NMC_all(1:6,:),2), mean(mu_NMC_all(1:6,:),2)+2*std(mu_NMC_all(1:6,:),0,2), mean(mu_NMC_all(1:6,:),2)-2*std(mu_NMC_all(1:6,:),0,2) ])



%%  RBF-BOED: multiple models


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 448 0.18; 180 0.21; 794 0.14];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(10,1);

mu_NMC_all = zeros(10,20);

N_M = 5;

Y2 = zeros(5000*N_M,61);

figure;

for j = 1:1

    Design = design_available(j,:);

    for k1 = 1:20

        for k = 1:6

            l_sum = 0;

            l = k*200-100;

            % l = 20;
            Y2_interp = [];
            % Y2 = zeros(5000*N_M,61);

            for i = 1:N_M

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                Model_prior = Model_design{i};

                model_true = Model_prior{1};

                mu_theta = Model_prior{3}.mu;

                sigma_theta = Model_prior{3}.sigma;

                data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
                    num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

                load(data_filename);

                % Y2((i-1)*5000+1:i*5000,:) = data{1}; 
                theta_target = data{2};


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                y2 = data{1};

                N_target = length(theta_target);

                k_std = 3.5;

                idx_all = 1:N_target;

                if i ==1 || i ==3

                    nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;


                     params_idx = [randi([1 length(y2)],l/4,1); (idx_all(boundary_idx ))'];

                      % params_idx = [(1:l/2)'; (idx_all(boundary_idx ))'];

                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                    [Nearest_Idx] = nearest_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],20);

                    [I]     =     RBF_PHS_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],Nearest_Idx,15,'1',3,1);

                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                     % y_interp = y2;


                elseif i==4

                    nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.25;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1;


                    nu_pos_idx_1 = theta_target(:,2)<1.5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                      % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                      params_idx  = unique( params_idx,'rows');

                      theta_params =  theta_target(params_idx,:);


                      % [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10],30);
                      % 
                      % [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10], Nearest_Idx,15,'1',3,1);



                      [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/15e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/5],[1.*(theta_params(:,1)/15e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/5],30);

                      [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/15e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/5],[1.*(theta_params(:,1)/15e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/5], Nearest_Idx,15,'1',3,1);


                      y_interp = I*y2(params_idx,:);
                      y_interp(params_idx,:) = y2(params_idx,:);

                      % y_interp = y2;

                      % filter_idx = true(size(y_interp,1),1);
                      %
                      % for j1 = 1:size(y_interp,1)
                      %
                      %     if max(y_interp(j1,:))>1+1e-6 || min(y_interp(j1,:))<min(min(y2(params_idx,:)))
                      %         filter_idx(j1) = false;
                      %     end
                    % 
                    % end
                    % y_interp =y_interp(filter_idx,:);

                     % [mu_NMC_j,~,] =  NMC_est(y_interp(:,1:61), [], sigma_w,'data');
                     % mean(mu_NMC_j)
                     % %
                     % 
                     % error = sqrt(mean(sigma_w(y2).*(y_interp-y2).^2,2)); figure; plot(error)
                     % 
                     % figure;
                     % 
                     % subplot(1,2,1)
                     % 
                     % target_idx = 4159;
                     % 
                     % scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,4),'k.'); hold on;
                     % scatter3(theta_target(target_idx,1),theta_target(target_idx,2),theta_target(target_idx,4),'ro'); hold on;
                     % scatter3(theta_params(Nearest_Idx(target_idx,1:15),1),theta_params(Nearest_Idx(target_idx,1:15),2),theta_params(Nearest_Idx(target_idx,1:15),4),'m*'); hold on;
                     % 
                     % 
                     % 
                     % axis square tight;
                     % 
                     % subplot(1,2,2)
                     % 
                     % plot(t,y_interp(target_idx,:),'r-'); hold on;
                     % 
                     % plot(t,y2(target_idx,:),'b--')
                     % 
                     % plot(t,y2(Nearest_Idx(target_idx,1:15),:),'k:')


                elseif i==5


                  nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.25;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1;


                    nu_pos_idx_1 = theta_target(:,2)<1.5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                  % [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10],30);
                  % 
                  %   [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10], Nearest_Idx,15,'1',3,1);
                  % 

                    [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/10e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/10e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6],30);

                    [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/10e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/10e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6], Nearest_Idx,15,'1',3,1);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);

                    % y_interp = y2;
                    % 
                    % [mu_NMC_j,~,] =  NMC_est(y_interp(:,1:61), [], sigma_w,'data');
                    %  mean(mu_NMC_j)
                    %  %
                    % 
                    %  error = sqrt(mean(sigma_w(y2).*(y_interp-y2).^2,2)); figure; plot(error)
                    % 
                    %  figure;
                    % 
                    %  subplot(1,2,1)
                    % 
                    %  target_idx = 3879;
                    % 
                    %  scatter3(theta_params(:,1),theta_params(:,2),theta_params(:,4),'k.'); hold on;
                    %  scatter3(theta_target(target_idx,1),theta_target(target_idx,2),theta_target(target_idx,4),'ro'); hold on;
                    %  scatter3(theta_params(Nearest_Idx(target_idx,1:15),1),theta_params(Nearest_Idx(target_idx,1:15),2),theta_params(Nearest_Idx(target_idx,1:15),4),'m*'); hold on;
                    % 
                    % 
                    % 
                    %  axis square tight;
                    % 
                    %  subplot(1,2,2)
                    % 
                    %  plot(t,y_interp(target_idx,:),'r-'); hold on;
                    % 
                    %  plot(t,y2(target_idx,:),'b--')
                    % 
                    %  plot(t,y2(Nearest_Idx(target_idx,1:15),:),'k:')
                    % 


                elseif i==2

                    % nan_idx = [173 279 371 1103 2060 2311 4406];
                    % 
                    %  y2(nan_idx,:) =  [];
                    % 
                    %  theta_target(nan_idx,:) =  [];

                    nu_pos_idx = theta_target(:,2)>5e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);


                    k_std= 3.55;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,3) -mu_theta(3)).^2/((sigma_theta(3,3))*k_std^2)>1;

                    params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))'];

                    % params_idx = [randi([1 N_target],l,1)];

                    % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                    % theta_target =  mvnrnd(mu_theta,sigma_theta,20000);


                    % [Nearest_Idx] = nearest_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],20);
                    %
                    % [I]     =     RBF_PHS_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],Nearest_Idx,20,'1',3,2);


                    [Nearest_Idx] = nearest_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],30);

                    [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],Nearest_Idx,27,'1',3,2);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);
                        % y_interp = y2;

                end

                % Y2_interp((i-1)*5000+1:i*5000,:) = y_interp;

                filter_idx = true(size(y_interp,1),1);

                for j1 = 1:size(y_interp,1)

                    for it = 1:61
                        if max(y_interp(j1,it))>(max(y2(params_idx,it)))+0.025 || min(y_interp(j1,it))<(min(y2(params_idx,it)))-0.025
                            filter_idx(j1) = false;
                        end
                    end

                end
                y_interp =y_interp(filter_idx,:);


                Y2_interp = [Y2_interp;y_interp];

                l_sum = l_sum+size( theta_params,1);

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [mu_NMC_j,~,] =  NMC_est(Y2_interp, [], sigma_w,'data');
            mu_NMC_all(k,k1)        =  mean(mu_NMC_j);

            mean(mu_NMC_j)

            plot(l_sum, mu_NMC_all(k,k1) ,'+'); hold on;

            axis square tight;

            xscale('log')

            xlim([10^2 5000*N_M]);

             ylim([4.7 6.6]);

             line([10 25000],[5.264 5.264],'Linestyle','--','Color','k')

            drawnow;

        end
    end


    clear data y2 theta_target


end


xlabel('We')
ylabel('Req')
zlabel('EIG')


csvwrite('IMR_5models_1_RBF.csv',[[400 1100 1700 2350 2950 3500]', mean(mu_NMC_all(1:6,1:20),2)-0.02, mean(mu_NMC_all(1:6,1:20),2)+2.5*std(mu_NMC_all(1:6,1:20),0,2)-0.02, mean(mu_NMC_all(1:6,1:20),2)-2.5*std(mu_NMC_all(1:6,1:20),0,2)-0.02 ])


%%


figure; 


load('EIG data/5models_EIG_1.mat');

plot(500*N_M*(1:10),mean(mu_NMC_all,2),'r.-'); hold on;

plot(500*N_M*(1:10),mean(mu_NMC_all,2)+2*std(mu_NMC_all,0,2),'r--'); hold on;

plot(500*N_M*(1:10),mean(mu_NMC_all,2)-2*std(mu_NMC_all,0,2),'r--'); hold on;

csvwrite('IMR_5models_2.csv',[ (500*N_M*(1:10))', mean(mu_NMC_all,2), mean(mu_NMC_all,2)+2*std(mu_NMC_all,0,2), mean(mu_NMC_all,2)-2*std(mu_NMC_all,0,2) ])

load('EIG data/5models_EIG_2_RBF.mat');

plot([500 1000 1500 2000 2500 3000],mean(mu_NMC_all(1:6,1:k1-1),2),'b.-'); hold on;

plot([500 1000 1500 2000 2500 3000],mean(mu_NMC_all(1:6,1:k1-1),2)+2*std(mu_NMC_all(1:6,1:k1-1),0,2),'b--'); hold on;

plot([500 1000 1500 2000 2500 3000],mean(mu_NMC_all(1:6,1:k1-1),2)-2*std(mu_NMC_all(1:6,1:k1-1),0,2),'b--'); hold on;


csvwrite('IMR_5models_0_RBF.csv',[[500 1000 1500 2000 2500 3000]', mean(mu_NMC_all(1:6,1:20),2), mean(mu_NMC_all(1:6,1:20),2)+3*std(mu_NMC_all(1:6,1:20),0,2), mean(mu_NMC_all(1:6,1:20),2)-3*std(mu_NMC_all(1:6,1:20),0,2) ])


   axis square tight;

        xlim([0 5000*N_M]);

    xscale('log')


% [500 800 1000 1200 1450]



%%  EIG convergence check: Multiple models


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 794 0.14; 448 0.18; 180 0.21; 204 0.25; 803 0.23];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(20,1);

mu_NMC_all = zeros(10,20);

N_M = 5;

Y2 = zeros(5000*N_M,61);

figure;

for j = 1:length(design_available)

    Design = design_available(j,:);

    Y2 =[];


    for i = 1:N_M

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

        filter_idx = true(size(y2,1),1);

        for j1 = 1:size(y2,1)

            for it = 1:61
                if max(y2(j1,it))>1.1 || min(y2(j1,it))<1e-6
                    filter_idx(j1) = false;
                end
            end

        end


         Y2 = [Y2;y2(filter_idx,:)];

    end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:20

           tic
                [mu_NMC_j,~,] =  NMC_est(Y2, [], sigma_w,'data');
                mu_NMC(k)        =  mean(mu_NMC_j);


            toc


            plot(k, mu_NMC(k) ,'.'); hold on;

            axis square tight;

            % xlim([0 5000*N_M]);

            % ylim([0.14 0.3]);

            drawnow;


        end

        clear data y2 theta_target


end


xlabel('We')
ylabel('Req')
zlabel('EIG')



csvwrite('IMR_1model_0_NHKV_RBF.csv',[ [33 134 233 333 433 532]', mean(mu_NMC_all(1:6,:),2), mean(mu_NMC_all(1:6,:),2)+2*std(mu_NMC_all(1:6,:),0,2), mean(mu_NMC_all(1:6,:),2)-2*std(mu_NMC_all(1:6,:),0,2) ])




% j = 4:  6.2817 +/- 2*0.0057
% j = 5:  6.1831 +/- 2*0.0075
% j = 6:  5.5289 +/- 2*0.0085
% j = 7:  5.5862 +/- 2*0.0071
% j = 8:  6.0854 +/- 2*0.0076

%%



design_available = [833 0.28; 475 0.26; 111 0.25; 794 0.14; 448 0.18; 180 0.21; 204 0.25; 803 0.23];

mu_NMC = [ normrnd(5.94,0.0057,1,1) normrnd(6.108,0.0075,1,1) normrnd(5.264,0.0085,1,1) ...
    normrnd(6.2817,0.0057,1,1) normrnd(6.1831,0.0075,1,1) normrnd(5.5289,0.0085,1,1) normrnd(5.5862,0.0071,1,1) normrnd(6.0854,0.0076,1,1)];


[Xgrid,Ygrid]    =  meshgrid(linspace(100,1000,51),...
    linspace(0.14,0.3,51));

  XYgrid              =  [Xgrid(:), Ygrid(:)];

  mdl = fitrgp(design_available,mu_NMC,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'));

  [y_pred,std]        =  predict(mdl,XYgrid);

  % mdl_interp = fitrgp(design_available,mu_NMC_interp(1:21),'KernelFunction','ardmatern52',...
  %           'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
  %           struct('AcquisitionFunctionName','expected-improvement-plus'));
  % 
  % [y_pred_interp,std_interp]        =  predict(mdl_interp,XYgrid);


figure; 

subplot(1,2,1)

surf(Xgrid,Ygrid,reshape(y_pred,51,51)); shading interp; hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC(:) ,'ro'); hold on;
   axis square tight;


    xlim([100 1000]);

    ylim([0.14 0.3]);

     zlim([5 6.5]);

    xlabel('We')
ylabel('Req')
zlabel('EIG')


% subplot(1,2,2)
% 
% surf(Xgrid,Ygrid,reshape(y_pred_interp,100,100)); shading interp; hold on;
% 
%     scatter3(design_available(:,1),design_available(:,2), mu_NMC_interp(1:21) ,'b*'); hold on;
%    axis square tight;
% 
% 
%     xlim([100 1000]);
% 
%     ylim([0.14 0.3]);
% 
%     zlim([0 3]);
% 
%     xlabel('We')
% ylabel('Req')
% zlabel('EIG')



csvwrite('EIG_field_5models.txt',[reshape([Ygrid; zeros(1,51)],[],1) reshape([Xgrid; zeros(1,51)],[],1) reshape([reshape(y_pred,51,51); zeros(1,51)],[],1)])


csvwrite('EIG_points_5models.csv',[design_available mu_NMC'])




%%  Bayesian optimization: single model


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 794 0.14; 448 0.18; 180 0.21;...
    204 0.25; 803 0.23; 250 0.15; 600 0.21; 1000 0.23; ...
      300 0.29; 350 0.14; 650 0.18; 900 0.17];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(length(design_available),1);


for j = 1:length(design_available)

    Design = design_available(j,:);


    for i = 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

       %  nu_pos_idx = theta_target(:,2)>0;
       % 
       % y2 = y2(nu_pos_idx,:);
       % theta_target = theta_target(nu_pos_idx,:);

        filter_idx = true(size(y2,1),1);

        for j1 = 1:size(y2,1)

            for it = 1:61
                if max(y2(j1,it))>1.1 || min(y2(j1,it))<1e-6
                    filter_idx(j1) = false;
                end
            end

        end


         y2 = y2(filter_idx,:);
           theta_target = theta_target(filter_idx,:);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         l = 500;

         N_target = length(theta_target);

         k_std = 3.5;

         idx_all = 1:N_target;

         if i ==1 || i ==3

             % nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;
             %
             % y2 = y2(nu_pos_idx,:);
             % theta_target = theta_target(nu_pos_idx,:);

             boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;


             params_idx = [randi([1 length(y2)],l,1); (idx_all(boundary_idx ))'];

             % params_idx = [(1:l/2)'; (idx_all(boundary_idx ))'];

             params_idx  = unique( params_idx,'rows');

             theta_params =  theta_target(params_idx,:);


             [Nearest_Idx] = nearest_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],20);

             [I]     =     RBF_PHS_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],Nearest_Idx,15,'1',3,1);

             y_interp = I*y2(params_idx,:);
             y_interp(params_idx,:) = y2(params_idx,:);


         elseif i==2

             % nan_idx = [173 279 371 1103 2060 2311 4406];
             %
             %  y2(nan_idx,:) =  [];
             %
             %  theta_target(nan_idx,:) =  [];

             nu_pos_idx = theta_target(:,2)>5e-2 & theta_target(:,1)>0;

             y2 = y2(nu_pos_idx,:);
             theta_target = theta_target(nu_pos_idx,:);


             k_std= 3.55;

             boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,3) -mu_theta(3)).^2/((sigma_theta(3,3))*k_std^2)>1;

             params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))'];

             % params_idx = [randi([1 N_target],l,1)];

             % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

             params_idx  = unique( params_idx,'rows');

             theta_params =  theta_target(params_idx,:);


             % theta_target =  mvnrnd(mu_theta,sigma_theta,20000);


             % [Nearest_Idx] = nearest_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],20);
             %
             % [I]     =     RBF_PHS_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],Nearest_Idx,20,'1',3,2);


             [Nearest_Idx] = nearest_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],30);

             [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],Nearest_Idx,27,'1',3,2);


             y_interp = I*y2(params_idx,:);
             y_interp(params_idx,:) = y2(params_idx,:);
             % y_interp = y2;

             filter_idx = true(size(y_interp,1),1);

             for j1 = 1:size(y_interp,1)

                 for it = 1:61
                     if max(y_interp(j1,it))>(max(y2(params_idx,it)))+0.005 || min(y_interp(j1,it))<(min(y2(params_idx,it)))-0.005
                         filter_idx(j1) = false;
                     end
                 end

             end
             y_interp =y_interp(filter_idx,:);

         end

     

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         [mu_NMC_j,~,] =  NMC_est(y2, [], sigma_w,'data');
         mu_NMC(j)        =  mean(mu_NMC_j);

         [mu_NMC_j,~,] =  NMC_est(y_interp, [], sigma_w,'data');
         mu_NMC_interp(j)        =  mean(mu_NMC_j);


    end

end



[Xgrid,Ygrid]    =  meshgrid(linspace(100,1000,51),...
    linspace(0.14,0.3,51));

  XYgrid              =  [Xgrid(:), Ygrid(:)];

  mdl = fitrgp(design_available,mu_NMC,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'));

  [y_pred,std]        =  predict(mdl,XYgrid);

  mdl_interp = fitrgp(design_available,mu_NMC_interp,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'));

  [y_pred_interp,std_interp]        =  predict(mdl_interp,XYgrid);


    xi                  =  -1;  % Exploration-exploitation parameter (High xi = more exploration, Low xi = more exploitation)
  deviation           =  y_pred - max(mu_NMC) - xi;
  EI                  =  (std ~= 0).*(std.*normpdf(deviation./std)+deviation.*normcdf(deviation./std));
  [~,EI_idx]          =  max(EI);
  x_next              =  XYgrid(EI_idx,:);


  deviation           =  y_pred_interp - max(mu_NMC_interp) - xi;
  EI                  =  (std_interp ~= 0).*(std_interp.*normpdf(deviation./std_interp)+deviation.*normcdf(deviation./std_interp));
  [~,EI_idx]          =  max(EI);
  x_next_interp       =  XYgrid(EI_idx,:);


figure; 

subplot(1,2,1)

surf(Xgrid,Ygrid,reshape(y_pred,51,51)); shading interp; hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC(:) ,'ro'); hold on;
   axis square tight;

       line([x_next(1) x_next(1)],[x_next(2),x_next(2)],[4 6.5])

    xlim([100 1000]);

    ylim([0.14 0.3]);

      zlim([4 6.5]);

    xlabel('We')
ylabel('Req')
zlabel('EIG')


subplot(1,2,2)

surf(Xgrid,Ygrid,reshape(y_pred_interp,51,51)); shading interp; hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC(:) ,'ro'); hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC_interp(:) ,'b*'); hold on;

        line([x_next_interp(1) x_next_interp(1)],[x_next_interp(2),x_next_interp(2)],[4 6.5])
   axis square tight;


    xlim([100 1000]);

    ylim([0.14 0.3]);

    zlim([4 6.5]);

    xlabel('We')
ylabel('Req')
zlabel('EIG')




csvwrite('EIG_field_qKV.txt',[reshape([Ygrid; zeros(1,51)],[],1) reshape([Xgrid; zeros(1,51)],[],1) reshape([reshape(y_pred,51,51); zeros(1,51)],[],1) ])


csvwrite('EIG_field_qKV_RBF.txt',[reshape([Ygrid; zeros(1,51)],[],1) reshape([Xgrid; zeros(1,51)],[],1) reshape([reshape(y_pred_interp,51,51); zeros(1,51)],[],1)])


csvwrite('EIG_points_qKV.csv',[design_available mu_NMC mu_NMC_interp'])




%%  Bayesian optimization: multi-model


t = linspace(0,3,61);

sigma_w = @(y) noise_fun(t,y);


design_available = [833 0.28; 475 0.26; 111 0.25; 794 0.14; 448 0.18; 180 0.21;...
    204 0.25; 803 0.23; 250 0.15; 600 0.21; 1000 0.23; ...
      300 0.29; 350 0.14; 650 0.18; 900 0.17];


norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

mu_NMC = zeros(length(design_available),1);



for j = 1:length(design_available)

    Design = design_available(j,:);

    Y2 = [];
Y2_interp = [];


    for i = 1:5

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
            num2str(Design(2),'%.2f') '/',model_true, '.mat']; % name of file containing R vs T data

        load(data_filename);

        y2 = data{1}; theta_target = data{2};

       %  nu_pos_idx = theta_target(:,2)>0;
       % 
       % y2 = y2(nu_pos_idx,:);
       % theta_target = theta_target(nu_pos_idx,:);

        filter_idx = true(size(y2,1),1);

        for j1 = 1:size(y2,1)

            for it = 1:61
                if max(y2(j1,it))>1.1 || min(y2(j1,it))<1e-6
                    filter_idx(j1) = false;
                end
            end

        end


         y2 = y2(filter_idx,:);
           theta_target = theta_target(filter_idx,:);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         l = 700;

         N_target = length(theta_target);

         k_std = 3.5;

         idx_all = 1:N_target;

         if i ==1 || i ==3

             % nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;
             %
             % y2 = y2(nu_pos_idx,:);
             % theta_target = theta_target(nu_pos_idx,:);

             boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)>1;


             params_idx = [randi([1 length(y2)],200,1); (idx_all(boundary_idx ))'];

             % params_idx = [(1:l/2)'; (idx_all(boundary_idx ))'];

             params_idx  = unique( params_idx,'rows');

             theta_params =  theta_target(params_idx,:);


             [Nearest_Idx] = nearest_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],20);

             [I]     =     RBF_PHS_interp([1.*(theta_target(:,1)/1.5e4) theta_target(:,2)],[1.*(theta_params(:,1)/1.5e4) theta_params(:,2)],Nearest_Idx,15,'1',3,1);

             y_interp = I*y2(params_idx,:);
             y_interp(params_idx,:) = y2(params_idx,:);


         elseif i==2

             % nan_idx = [173 279 371 1103 2060 2311 4406];
             %
             %  y2(nan_idx,:) =  [];
             %
             %  theta_target(nan_idx,:) =  [];

             nu_pos_idx = theta_target(:,2)>5e-2 & theta_target(:,1)>0;

             y2 = y2(nu_pos_idx,:);
             theta_target = theta_target(nu_pos_idx,:);


             k_std= 3.55;

             boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,3) -mu_theta(3)).^2/((sigma_theta(3,3))*k_std^2)>1;

             params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))'];

             % params_idx = [randi([1 N_target],l,1)];

             % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

             params_idx  = unique( params_idx,'rows');

             theta_params =  theta_target(params_idx,:);


             % theta_target =  mvnrnd(mu_theta,sigma_theta,20000);


             % [Nearest_Idx] = nearest_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],20);
             %
             % [I]     =     RBF_PHS_interp_3d([theta_target(:,[1 2]) tanh(theta_target(:,3)/2)],[theta_params(:,[1 2]) tanh(theta_params(:,3)/2)],Nearest_Idx,20,'1',3,2);


             [Nearest_Idx] = nearest_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],30);

             [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/3000 theta_target(:,2)*1 (theta_target(:,3)/3)],[theta_params(:,1)/3000 theta_params(:,2)*1 (theta_params(:,3)/3)],Nearest_Idx,27,'1',3,2);


             y_interp = I*y2(params_idx,:);
             y_interp(params_idx,:) = y2(params_idx,:);
             % y_interp = y2;

       

               elseif i==4

                    nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.25;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1;


                    nu_pos_idx_1 = theta_target(:,2)<1.5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                      % params_idx = [(1:l)'; (idx_all(boundary_idx ))'];

                      params_idx  = unique( params_idx,'rows');

                      theta_params =  theta_target(params_idx,:);


                      [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10],30);

                      [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10], Nearest_Idx,15,'1',3,1);



                      % [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/15e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/5],[1.*(theta_params(:,1)/15e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/5],30);
                      % 
                      % [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/15e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/5],[1.*(theta_params(:,1)/15e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/5], Nearest_Idx,15,'1',3,1);


                      y_interp = I*y2(params_idx,:);
                      y_interp(params_idx,:) = y2(params_idx,:);

                   

                elseif i==5


                  nu_pos_idx = theta_target(:,2)>1e-2 & theta_target(:,1)>0;

                    y2 = y2(nu_pos_idx,:);
                    theta_target = theta_target(nu_pos_idx,:);

                    idx_all = 1:length(theta_target);

                    k_std = 3.25;

                    boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1,1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2,2))*k_std^2)+(theta_target(:,4) -mu_theta(4)).^2/((sigma_theta(4,4))*k_std^2)>1;


                    nu_pos_idx_1 = theta_target(:,2)<1.5e-2;

                      params_idx = [randi([1 length(theta_target)],l,1); (idx_all(boundary_idx ))';idx_all(nu_pos_idx_1)' ];


                    params_idx  = unique( params_idx,'rows');

                    theta_params =  theta_target(params_idx,:);


                  [Nearest_Idx] = nearest_interp_3d( [theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(- (theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10],30);

                    [I]     =     RBF_PHS_interp_3d([theta_target(:,1)/2.5e4 theta_target(:,2) theta_target(:,2).*10.^(-(theta_target(:,4)-6)/2.75)/10],[theta_params(:,1)/2.5e4 theta_params(:,2) theta_params(:,2).*10.^(-(theta_params(:,4)-6)/2.75)/10], Nearest_Idx,15,'1',3,1);


                    % [Nearest_Idx] = nearest_interp_3d( [1.*(theta_target(:,1)/10e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1 -(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/10e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6],30);
                    % 
                    % [I]     =     RBF_PHS_interp_3d([1.*(theta_target(:,1)/10e4) theta_target(:,2)*1 (log10(theta_target(:,2))*1-(theta_target(:,4)-6)/1)/6],[1.*(theta_params(:,1)/10e4) theta_params(:,2)*1 (log10(theta_params(:,2))*1-(theta_params(:,4)-6)/1)/6], Nearest_Idx,15,'1',3,1);


                    y_interp = I*y2(params_idx,:);
                    y_interp(params_idx,:) = y2(params_idx,:);


         end

            filter_idx = true(size(y_interp,1),1);

             for j1 = 1:size(y_interp,1)

                 for it = 1:61
                     if max(y_interp(j1,it))>(max(y2(params_idx,it)))+0.000 || min(y_interp(j1,it))<(min(y2(params_idx,it)))-0.000
                         filter_idx(j1) = false;
                     end
                 end

             end
             y_interp =y_interp(filter_idx,:);



          Y2_interp = [Y2_interp;y_interp];

          Y2       = [Y2;y2];

     

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end


         [mu_NMC_j,~,] =  NMC_est(Y2, [], sigma_w,'data');
         mu_NMC(j)        =  mean(mu_NMC_j);

         [mu_NMC_j,~,] =  NMC_est(Y2_interp, [], sigma_w,'data');
         mu_NMC_interp(j)        =  mean(mu_NMC_j);


end



[Xgrid,Ygrid]    =  meshgrid(linspace(100,1000,51),...
    linspace(0.14,0.3,51));

  XYgrid              =  [Xgrid(:), Ygrid(:)];

  mdl = fitrgp(design_available,mu_NMC,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'));

  [y_pred,std]        =  predict(mdl,XYgrid);

  mdl_interp = fitrgp(design_available,mu_NMC_interp,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'));

  [y_pred_interp,std_interp]        =  predict(mdl_interp,XYgrid);



  xi                  =  0.01;  % Exploration-exploitation parameter (High xi = more exploration, Low xi = more exploitation)
  deviation           =  y_pred - max(mu_NMC) - xi;
  EI                  =  (std ~= 0).*(std.*normpdf(deviation./std)+deviation.*normcdf(deviation./std));
  [~,EI_idx]          =  max(EI);
  x_next              =  XYgrid(EI_idx,:);


  deviation           =  y_pred_interp - max(mu_NMC_interp) - xi;
  EI                  =  (std_interp ~= 0).*(std_interp.*normpdf(deviation./std_interp)+deviation.*normcdf(deviation./std_interp));
  [~,EI_idx]          =  max(EI);
  x_next_interp       =  XYgrid(EI_idx,:);



figure; 

subplot(1,2,1)

surf(Xgrid,Ygrid,reshape(y_pred,51,51)); shading interp; hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC(:) ,'ro'); hold on;

    line([x_next(1) x_next(1)],[x_next(2),x_next(2)],[4 6.5])
   axis square tight;


    xlim([100 1000]);

    ylim([0.14 0.3]);

      zlim([4 6.5]);

    xlabel('We')
ylabel('Req')
zlabel('EIG')


subplot(1,2,2)

surf(Xgrid,Ygrid,reshape(y_pred_interp,51,51)); shading interp; hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC(:) ,'ro'); hold on;

    scatter3(design_available(:,1),design_available(:,2), mu_NMC_interp(:) ,'b*'); hold on;

   line([x_next_interp(1) x_next_interp(1)],[x_next_interp(2),x_next_interp(2)],[4 6.5])
   axis square tight;


    xlim([100 1000]);

    ylim([0.14 0.3]);

    zlim([4 7.5]);

    xlabel('We')
ylabel('Req')
zlabel('EIG')




csvwrite('EIG_field_5models.txt',[reshape([Ygrid; zeros(1,51)],[],1) reshape([Xgrid; zeros(1,51)],[],1) reshape([reshape(y_pred,51,51); zeros(1,51)],[],1) ])


csvwrite('EIG_field_5models_RBF.txt',[reshape([Ygrid; zeros(1,51)],[],1) reshape([Xgrid; zeros(1,51)],[],1) reshape([reshape(y_pred_interp,51,51); zeros(1,51)],[],1)])


csvwrite('EIG_points_5models.csv',[design_available mu_NMC mu_NMC_interp'])






%%  process data


design_available = [833 0.28; 475 0.26; 111 0.25; 794 0.14; 448 0.18; 180 0.21;...
    204 0.25; 803 0.23; 250 0.15; 600 0.21; 1000 0.23; ...
      300 0.29; 350 0.14; 650 0.18; 900 0.17];



for j = 9

    Design = design_available(j,:);


    for i = 5

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Model_prior = Model_design{i};

        model_true = Model_prior{1};

        mu_theta = Model_prior{3}.mu;

        sigma_theta = Model_prior{3}.sigma;

        y2 = []; theta_target = [];

        for k1 = 1:2

            data_filename =  ['data_multimodels/', 'We',num2str(round(Design(1)),'%i'),'_Req', ...
                num2str(Design(2),'%.2f') '/',model_true, '_' num2str(k1) '.mat']; % name of file containing R vs T data

            load(data_filename);

            y2 = [y2; data{1}]; theta_target = [theta_target; data{2}];

        end

        data = {y2,theta_target};
        %

        save_path = ['We',num2str(round(Design(1)),'%i'),'_Req', ...
        num2str(Design(2),'%.2f')];

        data_filename =  ['data_multimodels/', save_path '/' model_true,'.mat']; % name of file containing R vs T data

        save(data_filename,'data','-v7.3')

    end

end

