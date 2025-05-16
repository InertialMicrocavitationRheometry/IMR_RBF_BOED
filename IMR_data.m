


addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')


%% Set up parallel enviroment

poolobj = parpool('Processes', 24);
fprintf('Number of workers: %g\n', poolobj.NumWorkers);


%%  Initialize the prior distribution

%%%%%%%%%%%%%%%%%%%%%%---Model 1: NHKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_prior_1.mu      =  [15090 0.209 0];
P_prior_1.sigma   =  ([4350 0 0; 0 0.18 0; 0 0 0]).^2;
model_prob_1      =  0.5;
Model_1_prior     =  {'NeoHook', 'G, mu, alpha', P_prior_1, model_prob_1};

% %%%%%%%%%%%%%%%%%%%%%%%---Model 2: qKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% P_prior_2.mu      =  [2770 0.286 0.28];
% P_prior_2.sigma   =  ([300 0 0; 0 0.186 0; 0 0 0.48]).^2;
% model_prob_2      =  0.5;
% Model_2_prior     =  {'fung', 'G_inf, mu, alpha', P_prior_2, model_prob_2};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theta_target =  mvnrnd(P_prior_1.mu,P_prior_1.sigma,5000);

count = 0;


while count<30
    %%  

   Design_opt = [unifrnd(100,1000), unifrnd(0.14,0.3)];

   disp(['Design #' num2str(count) ', Optimal design at We = ' num2str(round(Design_opt(1)),'%i') ', Req = ' num2str(Design_opt(2),'%.2f')])

    %% generate the data for the optimal design

    data_type     =  'sim'; % 'sim' or 'exp'
    data_set      =  'simulation_data';
    data_filepath =  'simulation_data/';
    data_filename =  ['sim_data', '_We',num2str(round(Design_opt(1)),'%i'),'_Req', ...
    num2str(Design_opt(2),'%.2f') '.mat']; % name of file containing R vs T data
    save_info     =  {data_type, data_set, data_filepath, data_filename};
    
    % default = 0.05, perform IMR simulations up to t = 4.

    [t,yth,Uth,dt] = IMR_simulations(theta_true,Design_opt,model_true,60,'auto',save_info);

    count = count+1;
end

