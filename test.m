

%% Clear everything 
close all; clear all; clc;
%%
N = 8^4;
addpath('./src/spectral/')

tic;

for j = 1:20
[t,R,U] = f_imrv2('stress',1,'radial',2,'nt',240,'mt',240); 

end

toc



%%


[t,R,U] = f_imrv2('stress',1,'radial',2,'mu',0.093,'g',8310);

figure; 
plot(t,R)