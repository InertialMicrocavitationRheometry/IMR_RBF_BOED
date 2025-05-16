addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')
addpath('functions/RBF')




%%    test2:3D


N = 2e4;

N_sample = 100;

% theta = [0.5 + 0.3.*randn(N,1) 0.3 + 0.7.*randn(N,1) 0.5 + 0.8.*randn(N,1)];

theta = [normrnd(0.5,0.3,N,1) normrnd(0.3,0.7,N,1) normrnd(0.5,0.8,N,1)];

theta_sample = theta(1:N_sample,:);

Nd = 11;

d = linspace(0,1,Nd);

mu_NMC = zeros(Nd,1);



mu1 = 0.1;
mu2 = -0.1;

sigma_1=0.05^2; sigma_2=0.05^2;

mu = [mu1;mu2];
sigma = cat(3,sigma_1,sigma_2);
p = ones(1,2)/2;
gm = gmdistribution(mu,sigma,p);

norm_fun = @(mu, x)       ( exp(-((x-mu-mu1).^2/(2*sigma_1)) )...
                             +exp(-((x-mu-mu2).^2/(2*sigma_2)) ));  % removed the constant here



CMAME2025_data = [1.83 2 2.137 2.115 2.105 2.106 2.119 2.139 2.169 2.202 2.250];

CMAME2025_voed_data = [ 
 1.808673469387755
 1.9744897959183674
 2.1096938775510203
2.0816326530612246
 2.066326530612245
 2.071428571428571
 2.079081632653061
 2.0969387755102042
 2.1224489795918364
2.158163265306123
2.196428571428572];




figure;

plot(d, CMAME2025_data,'b*'); hold on;

axis square tight;


ylim([1.7 2.3]);




for j = 1:Nd

     % for k = 1:20

         theta = [normrnd(0.5,0.3,N,1) normrnd(0.3,0.7,N,1) normrnd(0.5,0.8,N,1)];

         k_std = 4;

         theta_target = theta;

         mu_theta = [0.5 0.3 0.5];
         sigma_theta = [0.3^2 0.7^2 0.8^2];

         boundary_idx = (theta_target(:,1) -mu_theta(1)).^2/((sigma_theta(1))*k_std^2)+(theta_target(:,2) -mu_theta(2)).^2/((sigma_theta(2))*k_std^2)+(theta_target(:,3) -mu_theta(3)).^2/((sigma_theta(3))*k_std^2)>1;

         idx_all = 1:N;

         params_idx = [randi([1 N],N_sample,1); (idx_all(boundary_idx ))'];

         params_idx  = unique( params_idx,'rows');

         theta_sample =  theta_target(params_idx,:);

         [Nearest_Idx] = nearest_interp_3d(theta ,theta_sample,50);


         [I]     =     RBF_PHS_interp_3d(theta ,theta_sample,Nearest_Idx,10,'1',3,1);


         Y_sample =  theta_sample(:,1).^3*d(j)^2+theta_sample(:,2).*exp(-abs(0.2-d(j)))+sqrt(theta_sample(:,3).^2*d(j)*2);

        Y = I*Y_sample;

        Y0 =  theta(:,1).^3*d(j)^2+theta(:,2).*exp(-abs(0.2-d(j)))+sqrt(2*theta(:,3).^2*d(j)^1);

        error(j,:) = ((Y-Y0).^2);


        % [mu_NMC_j,~,] =  NMC_est_mixed(Y, [], [],'data');
        % mu_NMC(j,k)        =  mean(mu_NMC_j);

        % [mu_NMC_j0,~,] =  NMC_est_mixed(Y0, [], [],'data');
        % mu_NMC0(j,k)        =  mean(mu_NMC_j0);

        %  scatter(d(j),mu_NMC(j,k),'k.'); hold on;
        % 
        % scatter(d(j),mu_NMC0(j,k),'ro'); hold on;


        % scatter(d(j),error(j,:),'ro'); hold on;

        % 
        % ylim([1.7 2.3]);

        drawnow;

     % end

end


% xlabel('d'); ylabel('EIG')


csvwrite('toy_error_200_q1.csv',[200*ones(Nd-1,1), d(2:end)', log10([ mean(error(2:end,:),2)]), log10(mean(error(2:end,:),2)+2*std(error(2:end,:),0,2)),log10(mean(error(2:end,:),2)-2*std(error(2:end,:),0,2))]);



% 
% xlim([0 1]); 
% 
% ylim([0.5 2.5]);


figure; 

scatter3(theta(:,1),theta(:,2),theta(:,3),'k.');

axis equal;

hold on;

scatter3(theta_sample(:,1),theta_sample(:,2),theta_sample(:,3),'ro');


%% error analysis





%%


mu1 = 0.1;
mu2 = -0.1;

sigma_1=0.05^2; sigma_2=0.05^2;

mu = [mu1;mu2];
sigma = cat(3,sigma_1,sigma_2);
p = ones(1,2)/2;
gm = gmdistribution(mu,sigma,p);

norm_fun_log = @(mu, x)  log( 1/(2*sqrt(2*pi*sigma_1))*exp(-((x-mu-mu1).^2/(2*sigma_1)) )...
                             +1/(2*sqrt(2*pi*sigma_2))*exp(-((x-mu-mu2).^2/(2*sigma_2)) ));  % removed the constant here

% norm_fun_log = @(mu, Sigma, x)  +(-((x-mu)/(2*Sigma))*(x-mu)' );  % removed the constant here


    % y0 = 0.426;

 y0 = 1.235;

eps = 1e-6;

figure;

subplot(1,3,1)

S1_0 = exp(norm_fun_log(y0,Y0)).*exp(-((theta(:,1)-0.5).^2/(2*0.3^2)) );
S1   = exp(norm_fun_log(y0,Y)).*exp(-((theta(:,1)-0.5).^2/(2*0.3^2)) );

idx_0 = false(N,1);

idx_0(1) = true;

idx = false(N,1);

idx(1) = true;

        idx_true = 1;
for k = 2:N

    alpha = rand(1);

    if alpha<=(S1_0(k)/S1_0(idx_true))
        idx_0(k) = true;
                idx_true = k;
    end

    if alpha<=(S1(k)/S1(idx_true))
        idx(k) = true;
                idx_true = k;
    end
end

% scatter(theta(idx_0,1), S1_0(idx_0),'k.'); hold on;
% scatter(theta(idx,1), S1(idx),'ro'); hold on;

histogram(theta(idx_0,1),20,'Normalization','pdf'); hold on;
histogram(theta(idx,1),20,'Normalization','pdf');

sample1_0 = theta(idx_0,1);
sample1 = theta(idx,1);

xlim([-1 2])


subplot(1,3,2)

S2_0 = exp(norm_fun_log(y0,Y0)).*exp(-((theta(:,2)-0.3).^2/(2*0.7^2)) );
S2   = exp(norm_fun_log(y0,Y)).*exp(-((theta(:,2)-0.3).^2/(2*0.7^2)) );

idx_0 = false(N,1);

idx_0(1) = true;

idx = false(N,1);

idx(1) = true;

idx_true = 1;

for k = 2:N

    alpha = rand(1);

    if alpha<=(S2_0(k)/S2_0(idx_true))
        idx_0(k) = true;
        idx_true = k;
    end

    if alpha<=(S2(k)/S2(idx_true))
        idx(k) = true;
        idx_true = k;
    end
end

idx_0 = S2_0>eps;
idx = S2>eps;

histogram(theta(idx_0,2),20,'Normalization','pdf'); hold on;
histogram(theta(idx,2),20,'Normalization','pdf');

sample2_0 = theta(idx_0,2);
sample2 = theta(idx,2);

xlim([-2 2])


subplot(1,3,3)

S3_0 = exp(norm_fun_log(y0,Y0)).*exp(-((theta(:,3)-0.5).^2/(2*0.8^2)) );
S3   = exp(norm_fun_log(y0,Y)).*exp(-((theta(:,3)-0.5).^2/(2*0.8^2)) );

idx_0 = false(N,1);

idx_0(1) = true;

idx = false(N,1);

idx(1) = true;

        idx_true = 1;
for k = 2:N

    alpha = rand(1);

    if alpha<=(S3_0(k)/S3_0(idx_true))
        idx_0(k) = true;
                idx_true = k;
    end


    if alpha<=(S3(k)/S3(idx_true))
        idx(k) = true;
                idx_true = k;
    end
end


idx_0 = S2_0>eps;
idx = S2>eps;
% scatter(theta(idx_0,3), S3_0(idx_0),'k.'); hold on;
% scatter(theta(idx,3), S3(idx),'ro'); hold on;

histogram(theta(idx_0,3),20,'Normalization','pdf'); hold on;
histogram(theta(idx,3),20,'Normalization','pdf');

sample3_0 = theta(idx_0,3);
sample3 = theta(idx,3);

xlim([-3 3])



csvwrite('sample1_0_case3.csv',sample1_0);
csvwrite('sample1_case3.csv',sample1);
csvwrite('sample2_0_case3.csv',sample2_0);
csvwrite('sample2_case3.csv',sample2);
csvwrite('sample3_0_case3.csv',sample3_0);
csvwrite('sample3_case3.csv',sample3);

%%  plot data


NB = [50 100 200 500 2e4];

datasets = {'EIG data/toy_EIG_50.mat'; 'EIG data/toy_EIG_100.mat'; 'EIG data/toy_EIG_230.mat'; 'EIG data/toy_EIG_450.mat';'EIG data/toy_EIG_2e4.mat' };

Nd = 11;

d = linspace(0,1,Nd);

figure;

for k = 1:5

    data = load(datasets{k});

    mu_NMC_k = data.mu_NMC;

    % scatter3(d,mu_NMC_k,'.'); hold on;
    
    plot3(NB(k)*ones(Nd,1),d,mean(mu_NMC_k,2),'k.-'); hold on;

     plot3(NB(k)*ones(Nd,1),d,mean(mu_NMC_k,2)+2*std(mu_NMC_k,0,2),'k--'); hold on;
    plot3(NB(k)*ones(Nd,1),d,mean(mu_NMC_k,2)-2*std(mu_NMC_k,0,2),'k--'); hold on;

 xscale('log')

end


csvwrite('toy_EIG_50.csv',[NB(k)*ones(Nd,1), d', mean(mu_NMC_k,2), mean(mu_NMC_k,2)+2*std(mu_NMC_k,0,2),mean(mu_NMC_k,2)-2*std(mu_NMC_k,0,2)]);





