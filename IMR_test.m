


addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')


%% Set up an underlying model

model_true = 'neoHook';   

theta_true1 = [7310 0.093 0];

theta_true2 = [8310 0.063 0];

theta_true = (theta_true1+theta_true2)/2;

Design_opt = [480, 0.15];

nt = 81;

tic

[t,yth,Uth,dt] = IMR_simulations(theta_true,Design_opt,model_true,nt-1);   

[t,yth1,Uth1,dt] = IMR_simulations(theta_true1,Design_opt,model_true,nt-1); 

[t,yth2,Uth2,dt] = IMR_simulations(theta_true2,Design_opt,model_true,nt-1); 

toc

figure;

plot(t,yth,'r--'); hold on;

plot(t,yth1,'k.-'); hold on;

plot(t,yth2,'b-'); hold on;

plot(t,(yth1+yth2)/2,'mo-'); hold on;


%%  try Jacobian


model_true = 'neoHook';   
theta_true = [8310 0.093 0];

J = zeros(2,nt);

d_d1 = 1;
d_d2 = 2e-4;

% Design = [480, 0.15];

N_d1 = 20;
N_d2 = 20;

d1 = linspace(100,1000,N_d1);
d2 = linspace(0.14,0.3,N_d2);

J_all = zeros(N_d1,N_d2,2);

for id1 = 1:N_d1

    for id2 = 1:N_d2

        Design = [d1(id1),d2(id2)];

         [t,yth,Uth,~] = IMR_simulations(theta_true,[Design(1), Design(2)],model_true,80);

        [t,yth1,Uth1,~] = IMR_simulations(theta_true,[Design(1)-d_d1, Design(2)],model_true,80);

        [t,yth2,Uth2,~] = IMR_simulations(theta_true,[Design(1)+d_d1, Design(2)],model_true,80);

        J_all(id1,id2,1) = mean((yth2-yth1)/2/d_d1);

        [t,yth3,Uth3,~] = IMR_simulations(theta_true,[Design(1), Design(2)-d_d2],model_true,80);

        [t,yth4,Uth4,~] = IMR_simulations(theta_true,[Design(1), Design(2)+d_d2],model_true,80);

         J_all(id1,id2,2) = mean((yth4-yth3)/2/d_d2);

    end

end


[X,Y] = meshgrid(d1(1:end-4),d2(1:end-4));

figure; 

subplot(1,2,1)

pcolor(X,Y,squeeze(J_all(1:end-4,1:end-4,1))); axis square tight; shading interp;

colorbar;

subplot(1,2,2)

pcolor(X,Y,squeeze(J_all(1:end-4,1:end-4,2))); axis square tight; shading interp;

colorbar;



figure; 

plot(t,J_all(id1,id2,1),'r-','LineWidth',2); hold on;

plot(t,J_all(id1,id2,2),'b--','LineWidth',2); hold on;


figure; 

plot(t,yth,'r-'); hold on;

plot(t,yth4,'b--')

plot(t,yth3,'b--')



std = abs(yth-1)/25+(1-exp(-.1*t))/10;


figure; 

plot(t,yth+std,'r--'); hold on;

plot(t,yth,'b-'); hold on;

plot(t,yth-std,'r--');






