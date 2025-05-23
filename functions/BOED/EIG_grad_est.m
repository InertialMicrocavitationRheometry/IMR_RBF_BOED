function [grad_NMC_j,log_P_LKH] = EIG_grad_est(F,theta,sigma_w,varargin)
 
% nested Monte-Carlo estimator for the expected information gain

% Inputs: forward model, q = F(theta,d) + w
%         Prior distribution of theta, p(theta) or mu_theta and sigma_theta
%         design configuration, d
%         measurements y with dimension Nx * M;
%         distribution of the noise, p(w) or sigma_w
%         range for the model parameters theta_range = [theta_min thta_max];

%%  

%  Sample theta from the Prior distribution p(theta)
%  Calculate the results from the forward model


if nargin == 4
    F_type   =  varargin{1};
else
    F_type   =  'model';
end


if strcmpi(F_type,'model')
    q_model  =  F(theta);
elseif strcmpi(F_type,'data')
    q_model  =  F{1};
    J        =  F{2};
    dimension = size(J,2);
end


[N,Nx] = size(q_model);

norm_fun_log    = @(mu, Sigma, x)  log(1/sqrt(2*pi)^(length(Sigma))/det(Sigma))+(-((x-mu)/(2*Sigma))*(x-mu)' );  % removed the constant here

norm_grad_log = @(mu, Sigma, x,J)  log(1/sqrt(2*pi)^(length(Sigma))/det(Sigma))+(x-mu)/(1*Sigma)*transpose(J);  % removed the constant here


%%  draw samples using the likelihood function

y = zeros(N,Nx);
for n = 1:N
    y(n,:) = mvnrnd(q_model(n,:),sigma_w,1); 
end


%% Calculate the log Likelihood and its gradient

grad_log_P_LKH   = zeros(N,dimension);
grad_log_P_LKH_sum   = zeros(N,dimension);

log_P_LKH   = zeros(N,N);


for j1 = 1:N
     grad_log_P_LKH(j1,:) =   (norm_grad_log(q_model(j1,:), sigma_w,y(j1,:), J(j1,:)))';
    for j2 = 1:N

        log_P_LKH(j1,j2) =   (norm_fun_log(q_model(j2,:), sigma_w,y(j1,:)))';
        grad_log_P_LKH_sum(j1,:) = grad_log_P_LKH_sum(j1,:)+ (norm_grad_log(q_model(j1,:), sigma_w,y(j2,:), J(j1,:)))';
    end

end

P_EVD  = (sum(exp(log_P_LKH),2));

%% Calculate the gradient of EIG

grad_NMC_j   = ( grad_log_P_LKH-grad_log_P_LKH_sum./P_EVD);

end
