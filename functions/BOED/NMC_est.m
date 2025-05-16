function [mu_NMC_j,log_P_LKH] = NMC_est(F,theta,sigma_w,varargin)
 
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


if nargin >= 4
    F_type   =  varargin{1};
else
    F_type   =  'model';
end


if strcmpi(F_type,'model')
    q_model  =  F(theta);
elseif strcmpi(F_type,'data')
    q_model  =  F;
end

q_model0 = q_model;

if nargin >= 5

    q_model0   =  varargin{2};
end



[N,Nx] = size(q_model);

% norm_fun_log = @(mu, Sigma, x)  +(-((x-mu)/(2*Sigma))*(x-mu)' );  % removed the constant here

% norm_fun_log = @(mu, Sigma, x)  +(-(x-mu)./(sqrt(Sigma))*((x-mu)./sqrt(Sigma))' )/2;  % removed the constant here

norm_fun_log = @(mu, Sigma, x)  +-((x-mu)./(sqrt(Sigma)) ).^2/2;  % removed the constant here


%%  draw samples using the likelihood function

y = zeros(N,Nx);

SIGMA_W = zeros(N,Nx);


for n = 1:N

    if  isa(sigma_w,'function_handle')

        sigma_n = diag(sigma_w(q_model(n,:)));
    else

        sigma_n = sigma_w;
    end   

    SIGMA_W(n,:) = diag(sigma_n);

    y(n,:) = mvnrnd(q_model0(n,:),sigma_n,1);
end




%% Calculate the likelihood


log_P_LKH   = zeros(N,N);

% for j1 = 1:N
%     for j2 = 1:N
%         log_P_LKH(j1,j2) =   (norm_fun_log(q_model(j2,:), (SIGMA_W(j2,:)),y(j1,:)))';
%     end
% end


    for j2 = 1:N
        % % 
        %  size((norm_fun_log(q_model(j2,:), (SIGMA_W(j2,:)),y))');

        log_P_LKH(:,j2) =   (sum(norm_fun_log(q_model(j2,:), (SIGMA_W(j2,:)),y),2))';
    end


log_P_EVD  = log(mean(exp(log_P_LKH),2));

%% Calculate the EIG

mu_NMC_j   = diag(log_P_LKH) -log_P_EVD;

end
