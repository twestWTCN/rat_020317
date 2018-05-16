function [xq yq R2] = vmfit(x,y,xn,plotop)
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
fit = @(b,x) (exp(b(1).*cos(x-b(2))))./(2*pi*b(3)*besselj(0,b(1);    % Function to fit

fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
options = optimset('MaxFunEvals',1e6);
[s] = fminsearchcon(fcn, [yr;  per;  -1;  ym],[-300; 0.75*pi; -1; -500],[300;2*pi; pi; 500],[],[],[],options)       ;                % Minimise Least-Squares
yfit = fit(s,x);
xq = linspace(min(x),max(x),xn);
yq = fit(s,xq);

if plotop == 1
    plot(x,y,'b',  xq,yq, 'r')
    grid
end
SSres = (y-yfit).^2;
SStot = (y-mean(y)).^2; % total sum of squares (variance)
R2 = sum(SSres)/sum(SStot);


