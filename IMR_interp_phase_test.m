



y_avail =y2(params_idx,:);


locs_j = zeros(size(y_avail,1),2);

y_phase_sep = zeros(size(y_avail,1),120+4);


for j2 = 1:size(y_avail,1)

[~,locs] = findpeaks(1-y_avail(j2,:));


locs_j(j2,1:2) = locs(1:2);


 y_phase_sep(j2,1:40)      =  interp1(t(1:locs(1)),y_avail(j2,1:locs(1)), linspace(t(1),t(locs(1)),40), 'makima' );

 y_phase_sep(j2,41:80)     =  interp1(t(locs(1)+1:locs(2)),y_avail(j2,locs(1)+1:locs(2)), linspace(t(locs(1)+1),t(locs(2)),40), 'makima' );

 y_phase_sep(j2,81:120)     =  interp1(t(locs(2)+1:end),y_avail(j2,locs(2)+1:end), linspace(t(locs(2)+1),t(end),40), 'makima' );

  y_phase_sep(j2,121:124)  = t([locs(1) locs(1)+1 locs(2) locs(2)+1]); 

end


y_interp0 = I*y_phase_sep;


for j1 = 1:size(y_interp0,1)

    y_interp_2(j1,:) =  interp1([linspace(t(1),y_interp0(j1,121),40) linspace(y_interp0(j1,122),y_interp0(j1,123),40) linspace(y_interp0(j1,124),t(end),40)],y_interp0(j1,1:120),t, 'makima');

end



norm_fun_log = @(mu, Sigma, x)  +(((x-mu)/(1*Sigma))*(x-mu)' );  % removed the constant here

 for k = 1:N_target
            Sigma = diag(sigma_w(y2(k,:)));
            % Sigma = eye(61,61);
            error(k) = (norm_fun_log(y2(k,:), Sigma,y_interp(k,:)))/1^2;
            error_2(k) = (norm_fun_log(y2(k,:), Sigma,y_interp_2(k,:)))/1^2;

 end


figure; 
plot(error,'b.-'); hold on;
plot(error_2,'r.-');



target_idx = 1184;


[t0,y_target,~,~] = IMR_simulations(theta_target(target_idx,:),Design,model_true,120);





sigma_w = @(y) noise_fun(t0,y);

figure; 

subplot(2,1,1)

plot(t0,y_target,'ro-'); hold on;

plot(t,y_interp(target_idx,:),'b*-'); hold on;

% plot(t,y_interp_2(target_idx,:),'m.-'); hold on;

  plot(t,y_interp_2(target_idx,:),'m+-'); hold on;

plot(t0,y_target+sqrt(sigma_w(y_target))*2,'r--');
plot(t0,y_target-sqrt(sigma_w(y_target))*2,'r--');

 plot(t,y2(Nearest_Idx(target_idx,1:20),:),'k:'); hold on;

  % plot( y_phase_sep(Nearest_Idx(target_idx,1:20),1:60)','k:'); hold on;

ylim([0 1])

legend('target','interpolation','interpolation (phase-locked)')

xlabel('t')

ylabel('R')


subplot(2,1,2)

scatter(theta_params(:,1),theta_params(:,2),'b.'); axis square tight; hold on; 

scatter(theta_target(target_idx,1), theta_target(target_idx,2),'r*');

scatter(theta_params(Nearest_Idx(target_idx,1:20),1),theta_params(Nearest_Idx(target_idx,1:20),2),'ko');

xlabel('G');
ylabel('\mu')

legend('samples','target','neighbors')


%%

figure;


plot( y_phase_sep(Nearest_Idx(target_idx,1:20),1:120)','ko-'); hold on;

% plot(y_interp0(target_idx,1:120),'r-');

plot(I(target_idx,:)*y_phase_sep(:,1:120),'r-')

