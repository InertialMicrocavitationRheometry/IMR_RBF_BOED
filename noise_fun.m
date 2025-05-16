function [sigma] = noise_fun(t,y)


sigma = (abs(y-1)/20+(1-exp(-t/10))/30).^2;

sigma(1) = 1e-10;

end