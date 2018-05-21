function [dx] = sgolaydiff(x,fsamp,npoly,framez)
if nargin<3
    npoly = 5;
end
if nargin<4
    framez = 25;
end

dt = 1/fsamp;
[b,g] = sgolay(npoly,framez);
dx = zeros(length(x),4);

p =1; % 1st order derivative
dx = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
