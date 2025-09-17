%% Settings

discreteZn = true; % Plot the purple triangles? Requires Chebfun.
output_tikz = true; % Output tikz files for the figures? Requires matlab2tikz.
makemultiplefigures = true; % Make 4 separate, new figures?

%% Helper functions

function A = sample_kernel(K, xpts, ypts)
A = zeros(length(xpts),length(ypts));
for i = 1:length(xpts), for j = 1:length(ypts)
    A(i,j) = K(xpts(i),ypts(j));
end, end
end
% Returns the weights for the first barycentric formula:
function w = firstBaryWeights(roots, poles)
n = length(roots);
w = ones(n,1);
for j = 1:n
    vp = roots(j) - poles;
    vr = roots(j) - roots;
    vr(j) = 1;
    w(j) = exp(sum(log(abs(vp)))-sum(log(abs(vr))));
    w(j) = w(j) * (prod(sign(vp))/prod(sign(vr)));
end
end
% Evaluates the first barycentric interpolation formula:
function Kn = firstBary(K,x,y,pts,pls,wts)
if isempty(pts)
    Kn = 0;
elseif any(pts==y)
    Kn = K(x,y);
else
    phi = exp(sum(log(abs(y-pts)))-sum(log(abs(y-pls))));
    phi = phi * (prod(sign(y-pts)) / prod(sign(y-pls)));
    Kn = phi * sum(K(x,pts).*wts./(y-pts));
end
end
% Calculates the error between the sampled matrix and the low rank
% approximation:
function err = interp_error(K,xpts,ypts,pts,pls,wts)
E = zeros(length(xpts),length(ypts));
for i = 1:length(xpts), for j = 1:length(ypts)
    E(i,j) = K(xpts(i),ypts(j)) - firstBary(K,xpts(i),ypts(j),pts,pls,wts);
end, end
err = norm(E);
end

%% Computation for Figure 1

N = 100;
s = 1/2;
t = 1/2;
K = @(x,y) beta(x+y+s,t);
xpts = (0:N)';
ypts = (0:N)';
A = sample_kernel(K, xpts, ypts);

svals = svd(A);

a = -Inf; b = -s; c = 0; d = N;

function out = w(k,t)
    if k + 1 - t <= 0 % gammaln needs positive input
        out = gamma(k+1-t)./(gamma(k+1));
    else
        out = exp(gammaln(k+1-t) - gammaln(k+1));
    end
    out = abs(out/gamma(1-t));
end

Kpnorm = norm(w(0:N,t)/w(N,t));

Cnorm = 0;
for y = 0:N
    tmp = sum(w(0:(2*N),t)./(y+s+(0:(2*N))));
    k = 2*N+1;
    toadd = w(k,t)/(y+k+s);
    while toadd/tmp > 1e-4
        tmp = tmp + toadd;
        k = k + 1;
        toadd = w(k,t)/(y+k+s);
    end
    Cnorm = Cnorm + tmp^2;
end
Cnorm = sqrt(Cnorm);

r = 40;
Z = ZolotarevNumber(a,b,c,d,0:r);
chebinterperr = zeros(r+1,1);
zolinterperr = zeros(r+1,1);
for n = 0:r
    pts = chebpts(n, [c,d], 1);
    wts = firstBaryWeights(pts,[]);
    chebinterperr(n+1) = interp_error(K,xpts,ypts,pts,[],wts);
    [pts,pls] = ellipjZolotarev(a,b,c,d,n);
    wts = firstBaryWeights(pts,pls);
    zolinterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
end
chebinterperr = chebinterperr/svals(1);
zolinterperr = zolinterperr/svals(1);

R = 1 + (2/(d-c)) * ( c - b + sqrt((d-b)*(c-b)));
% F is real, so use Remark 5.4:
BernsteinEst = Cnorm*Kpnorm*(2*R.^(-(0:r)))./(1+R.^(-2*(0:r))) / svals(1);

BTinds = floor((0:r)/2) - mod(N,2);
BTEst = 16*exp((-BTinds)*pi^2 / (2*log(8*floor(N/2)/pi)));

WebbEst = Z.*Cnorm*Kpnorm/svals(1);

if discreteZn
    aaaZ = zeros(r+1,1);
    aaainterperr = zeros(r+1,1);
    E = ypts;
    % Use the elements of b, b-1, b-2, b-3,...
    % that are closest to well-distributed points
    [~,F] = ellipjZolotarev(a,b,c,d,1000);
    F = unique(b + round(F-b));
    F = [-Inf; F];
    toInf = d+1;  % Mobius maps d+1 to Inf.
    TE = 1./(toInf-E);
    TF = 1./(toInf-F);
    for n = 12:22
        [pts, pls, sigma] = aaaZolotarev(TE,TF,n);
        aaaZ(n+1) = sigma;
        pts = toInf - 1./pts; % invert Mobius transformation
        pls = toInf - 1./pls;
        wts = firstBaryWeights(pts,pls);
        aaainterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
    end
    aaaWebbEst = aaaZ.*Cnorm*Kpnorm/svals(1);
    aaainterperr = aaainterperr/svals(1);
end

svals = svals/svals(1);

%% Plot Figure 1

if makemultiplefigures, figure, end

zolinterperr(27:end) = 0;
aaainterperr(21:end) = 0;

colours = orderedcolors("gem");
semilogy(0:r, BernsteinEst, 'color', colours(1,:))
hold on
semilogy(0:r, chebinterperr, '.', 'color', colours(1,:),'MarkerSize',5)
semilogy(0:r, BTEst, 'color', colours(2,:))
semilogy(0:r, WebbEst, 'color', colours(3,:))
semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
if discreteZn
    semilogy(0:r, aaaWebbEst, 'color', colours(4,:))
    semilogy(0:r, aaainterperr, '^', 'color', colours(4,:))
end
semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
hold off;

xlabel('$n$ (rank)', 'Interpreter', 'latex')
ylabel('$\|A - A_n \|_2 / \|A\|_2$', 'Interpreter', 'latex')
if discreteZn
    legend("Little--Reade \cite{little1984eigenvalues}", 'Chebyshev interpolant', "Beckermann--Townsend \cite{beckermann2019bounds}", 'Theorem \ref{thm:main}', 'Zolotarev interpolant', "Theorem \ref{thm:main} (discrete)", 'Discrete Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
else
    legend("Little--Reade \cite{little1984eigenvalues}", 'Chebyshev interpolant', "Beckermann--Townsend \cite{beckermann2019bounds}", 'Theorem \ref{thm:main}', 'Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
end
axis([0,r,1.25e-14, 1e1])
set(findall(gcf,'Type','line'),'LineWidth',0.8);

cleanfigure
if output_tikz
    matlab2tikz('figure1.tex', 'height', '3.5cm', 'width', '7.5cm');
end

%% Computation and plot for Figure 3 left
% Cauchy Matrices

N = 100;
K = @(x,y) 1./(y+x);

D = [1, 70];
E = [2,100];
xpts = linspace(D(1),D(2),N)';
ypts = linspace(E(1),E(2),N)';
A = sample_kernel(K, xpts, ypts);

svals = svd(A);

a = -xpts(end); b = -xpts(1); c = ypts(1); d = ypts(end);

Kpnorm = 1;
Cnorm = svals(1);

r = 25;
Z = ZolotarevNumber(a,b,c,d,0:r);
chebinterperr = zeros(r+1,1);
zolinterperr = zeros(r+1,1);
for n = 0:r
    pts = chebpts(n, [c,d], 1);
    wts = firstBaryWeights(pts,[]);
    chebinterperr(n+1) = interp_error(K,xpts,ypts,pts,[],wts);
    [pts,pls] = ellipjZolotarev(a,b,c,d,n);
    wts = firstBaryWeights(pts,pls);
    zolinterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
end
chebinterperr = chebinterperr/svals(1);
zolinterperr = zolinterperr/svals(1);

R = 1 + (2/(d-c)) * ( c - b + sqrt((d-b)*(c-b)));
% F is real, so use Remark 5.4
BernsteinEst = Kpnorm*Cnorm*(2*R.^(-(0:r)))./(1+R.^(-2*(0:r)))/svals(1);

WebbEst = Z.*Kpnorm*Cnorm/svals(1);

svals = svals/svals(1);

% Figure 3 left

if makemultiplefigures, figure, end

nsp = 13;

subplot(1,2,1)

zolinterperr(20:end) = 0;

colours = orderedcolors("gem");

semilogy(0:r, BernsteinEst, 'color', colours(1,:))
hold on
semilogy(0:r, chebinterperr, '.', 'color', colours(1,:),'MarkerSize',5)
semilogy(0:r, WebbEst, 'color', colours(3,:))
semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
hold off;

xlabel('$n$ (rank)', 'Interpreter', 'latex')
ylabel('$\|A - A_n \|_2 / \|A\|_2$', 'Interpreter', 'latex')
axis([0,r,1e-14, 1e1])


% Computation for Figure 3 right

N = 50;
K = @(x,y) 1./(x+y);

D = [1, 70, 1, 199];
E = [2, 100];
F = [-(D(2)+D(4)),-(D(1)+D(3))];
wpts = linspace(D(1),D(2),N)';
xpts = linspace(D(3),D(4),N)';
wxpts = wpts + xpts';
wxpts = wxpts(:);
ypts = linspace(E(1),E(2),N)';
A = sample_kernel(K, wxpts, ypts);

svals = svd(A);

a = F(1); b = F(2); c = E(1); d = E(2);

Kpnorm = 1;
Cnorm = svals(1);

r = 25;
Z = ZolotarevNumber(a,b,c,d,0:r);
chebinterperr = zeros(r+1,1);
zolinterperr = zeros(r+1,1);
for n = 0:r
    pts = chebpts(n, [c,d], 1);
    wts = firstBaryWeights(pts,[]);
    chebinterperr(n+1) = interp_error(K,wxpts,ypts,pts,[],wts);
    [pts,pls] = ellipjZolotarev(a,b,c,d,n);
    wts = firstBaryWeights(pts,pls);
    zolinterperr(n+1) = interp_error(K,wxpts,ypts,pts,pls,wts);
end
chebinterperr = chebinterperr/svals(1);
zolinterperr = zolinterperr/svals(1);

R = 1 + (2/(d-c)) * ( c - b + sqrt((d-b)*(c-b)));
% F is real, so use Remark 5.4
BernsteinEst = Kpnorm*Cnorm*(2*R.^(-(0:r)))./(1+R.^(-2*(0:r)))/svals(1);

WebbEst = Z.*Kpnorm*Cnorm/svals(1);

if discreteZn
    aaaZ = zeros(r+1,1);
    aaainterperr = zeros(r+1,1);
    E = ypts;
    F = -wxpts;
    for n = 13:19
        [pts, pls, sigma] = aaaZolotarev(E,F,n);
        aaaZ(n+1) = sigma;
        wts = firstBaryWeights(pts,pls);
        aaainterperr(n+1) = interp_error(K,wxpts,ypts,pts,pls,wts);
    end
    aaaWebbEst = aaaZ.*Kpnorm*Cnorm/svals(1);
    aaainterperr = aaainterperr/svals(1);
end

svals = svals/svals(1);


subplot(1,2,2)


aaainterperr(18:end) = 0;
zolinterperr(21:end) = 0;

colours = orderedcolors("gem");

semilogy(0:r, BernsteinEst, 'color', colours(1,:))
hold on
semilogy(0:r, chebinterperr, '.', 'color', colours(1,:),'MarkerSize',5)
semilogy(0:r, WebbEst, 'color', colours(3,:))
semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
if discreteZn
    semilogy(0:r, aaaWebbEst, 'color', colours(4,:))
    semilogy(0:r, aaainterperr, '^', 'color', colours(4,:))
end
semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
hold off;

xlabel('$n$ (rank)', 'Interpreter', 'latex')

axis([0,r,1e-14, 1e1])

if discreteZn
    legend("Proposition \ref{prop:LR}", 'Chebyshev interpolant', 'Theorem \ref{thm:main}', 'Zolotarev interpolant', "Theorem \ref{thm:main} (discrete)", 'Discrete Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
else
    legend("Proposition \ref{prop:LR}", 'Chebyshev interpolant', 'Theorem \ref{thm:main}', 'Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
end

set(gca,'yticklabel',[])

set(findall(gcf,'Type','line'),'LineWidth',0.8);

cleanfigure
if output_tikz
    matlab2tikz('figure3.tex', 'height', '3.5cm', 'width', '10cm');
end




%% Computation for Figure 4
% LogKernel

N = 100;
K = @(x,y) log(x+y);
c = 1; d = N;
xpts = c + [0;rand(N-2,1);1]*(d-c);
ypts = c + [0;rand(N-2,1);1]*(d-c);
A = sample_kernel(K, xpts, ypts);

svals = svd(A);

c = ypts(1); d = ypts(end);
a = -Inf; b = -c;

Kpnorm = sqrt(N);
Cnorm = 0.5*log((d+c)/(2*c))*sqrt(N); % the pre-factor in Lemma 3.5.

r = 40;

Z = ZolotarevNumber(a,b,c,d,0:r);
chebinterperr = zeros(r+1,1);
zolinterperr = zeros(r+1,1);
zolinterperr2 = zeros(r+1,1);
for n = 0:r-1
    pts = chebpts(n, [c,d], 1);
    pts = [pts; b + sqrt((d-b)*(c-b))];
    wts = firstBaryWeights(pts,[]);
    chebinterperr(n+2) = interp_error(K,xpts,ypts,pts,[],wts);
    [pts,pls,sigma] = ellipjZolotarev(a,b,c,d,n);
    pts = [pts; b + sqrt((d-b)*(c-b))+mod(n,2)/2]; % avoid repeated point
    pts = sort(pts);
    wts = firstBaryWeights(pts,pls);
    Z(n+2) = sigma;
    zolinterperr(n+2) = interp_error(K,xpts,ypts,pts,pls,wts);
    [pts,pls] = ellipjZolotarev(a,b,c,d,n);
    wts = firstBaryWeights(pts,pls);
    zolinterperr2(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
end
chebinterperr = chebinterperr/svals(1);
zolinterperr = zolinterperr/svals(1);
zolinterperr2 = zolinterperr2/svals(1);

R = 1 + (2/(d-c)) * ( c - b + sqrt((d-b)*(c-b)));
% F is real, so use Remark 5.4
BernsteinEst = Kpnorm*Cnorm*(2*R.^(-(0:(r-1))))./(1+R.^(-2*(0:(r-1))))/svals(1);
BernsteinEst = [0;BernsteinEst.'];

WebbEst = Z.*Kpnorm*Cnorm/svals(1);

svals = svals/svals(1);

%% Plot Figure 4

if makemultiplefigures, figure, end

zolinterperr(22:end) = 0;
zolinterperr2(25:end) = 0;

colours = orderedcolors("gem");
semilogy(0:r, BernsteinEst, 'color', colours(1,:))
hold on
semilogy(0:r, chebinterperr, '.', 'color', colours(1,:),'MarkerSize',5)
semilogy(0:r, WebbEst, 'color', colours(3,:))
semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
semilogy(0:r, zolinterperr2, '+', 'color', colours(6,:))
semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
hold off;

xlabel('$n$ (rank)', 'Interpreter', 'latex')
ylabel('$\|A - A_n \|_2 / \|A\|_2$', 'Interpreter', 'latex')

legend("Proposition \ref{prop:LR}", 'Chebyshev interpolant', "Equation \eqref{eqn:logCauchybound}", 'Zolotarev interpolant', 'Suboptimal Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');

axis([0,r,1e-14, 1e1])
set(findall(gcf,'Type','line'),'LineWidth',0.8);

cleanfigure
if output_tikz
    matlab2tikz('figure4.tex', 'height', '3.5cm', 'width', '7.5cm');
end


%% Computation for Figure 5
% Twisted Hankel transform matrix

function rt = besselroots(n)
    rt = zeros(size(n));
    for k = 1:length(n)
        rt(k) = besselroot(n(k));
    end
end
function rt = besselroot(n)
    BESSELJ0_ROOTS = [2.4048255576957728, 5.5200781102863106, 8.6537279129110122, 11.791534439014281, ...
    14.930917708487785, 18.071063967910922, 21.211636629879258, 24.352471530749302, ...
    27.493479132040254, 30.634606468431975, 33.775820213573568, 36.917098353664044, ...
    40.058425764628239, 43.199791713176730, 46.341188371661814, 49.482609897397817, ...
    52.624051841114996, 55.765510755019979, 58.906983926080942, 62.048469190227170];
    if n <= 20
        rt = BESSELJ0_ROOTS(n);
    else
        b = pi*(n-1/4);
        a1 = 1/8;
        a3 = -31/384;
        a5 = 3779/15360;
        a7 = -6277237/3440640;
        a9 = 2092163573/82575360;
        a11 = -8249725736393/14533263360;
        a13 = 423748443625564327/22671890841600;
        rt = b + ((((((a13/b^2 + a11)/b^2 + a9)/b^2 + a7)/b^2 + a5)/b^2 + a3)/b^2 + a1)/b;
    end
end

N = 100;
K = @(x,y) besselh(0,x.*y).*exp(-1i*x.*y);
xpts = besselroots((1:N)')/besselroot(N+1);
ypts = besselroots(1:N)';
A = sample_kernel(K, xpts, ypts);

svals = svd(A);

a = -Inf; b = 0; c = ypts(1); d = ypts(end);

Kpnorm = (2/pi) * sqrt(N);
Cnorm = 0.5*log((d+c)/(2*c))*sqrt(N);

r = 40;

Z = ZolotarevNumber(a,b,c,d,0:r);
chebinterperr = zeros(r+1,1);
zolinterperr = zeros(r+1,1);
zolinterperr2 = zeros(r+1,1);
for n = 0:r-1
    pts = chebpts(n, [c,d], 1);
    pts = [pts; b + sqrt((d-b)*(c-b))];
    wts = firstBaryWeights(pts,[]);
    chebinterperr(n+2) = interp_error(K,xpts,ypts,pts,[],wts);
    [pts,pls,sigma] = ellipjZolotarev(a,b,c,d,n);
    pts = [pts; b + sqrt((d-b)*(c-b))+mod(n,2)/2]; % avoid repeated point
    pts = sort(pts);
    wts = firstBaryWeights(pts,pls);
    Z(n+2) = sigma;
    zolinterperr(n+2) = interp_error(K,xpts,ypts,pts,pls,wts);
end
chebinterperr = chebinterperr/svals(1);
zolinterperr = zolinterperr/svals(1);

R = 1 + (2/(d-c)) * ( c - b + sqrt((d-b)*(c-b)));
% F is real, so use Remark 5.4
BernsteinEst = Kpnorm*Cnorm*(2*R.^(-(0:(r-1))))./(1+R.^(-2*(0:(r-1))))/svals(1);
BernsteinEst = [0;BernsteinEst.'];

WebbEst = Z.*Kpnorm*Cnorm/svals(1);

if discreteZn
    aaaZ = zeros(r+1,1);
    aaainterperr = zeros(r+1,1);
    E = ypts;
    toInf = d+1;  % Mobius maps d+1 to Inf.
    TE = 1./(toInf-E);
    TF = linspace(0,1/toInf,2000)';
    for n = 12:22
        [pts, pls, sigma] = aaaZolotarev(TE,TF,n);
        aaaZ(n+1) = sigma;
        pts = toInf - 1./pts; % invert Mobius transformation
        pls = toInf - 1./pls;
        wts = firstBaryWeights(pts,pls);
        aaainterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
    end
    aaaWebbEst = aaaZ.*Cnorm*Kpnorm/svals(1);
    aaainterperr = aaainterperr/svals(1);
end

svals = svals/svals(1);

%% Plot Figure 5

zolinterperr(25:end) = 0;
aaainterperr(21:end) = 0;

if makemultiplefigures, figure, end

colours = orderedcolors("gem");
semilogy(0:r, BernsteinEst, 'color', colours(1,:))
hold on
semilogy(0:r, chebinterperr, '.', 'color', colours(1,:),'MarkerSize',5)
semilogy(0:r, WebbEst, 'color', colours(3,:))
semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
if discreteZn
    semilogy(0:r, aaainterperr, '^', 'color', colours(4,:))
end
semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
hold off;

xlabel('$n$ (rank)', 'Interpreter', 'latex')
ylabel('$\|A - A_n \|_2 / \|A\|_2$', 'Interpreter', 'latex')
if discreteZn
    legend("Proposition \ref{prop:LR}", 'Chebyshev interpolant', "Equation \eqref{eqn:Hankelbound}", 'Zolotarev interpolant', 'Semidiscrete Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
else
    legend("Proposition \ref{prop:LR}", 'Chebyshev interpolant', "Equation \eqref{eqn:Hankelbound}", 'Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
end
axis([0,r,1.25e-14, 1e1])
set(findall(gcf,'Type','line'),'LineWidth',0.8);

cleanfigure
if output_tikz
    matlab2tikz('figure5.tex', 'height', '3.5cm', 'width', '7.5cm');
end




%% Not included in the paper because we exactly recover the best approximation
% 
% % circle with centre a and radius b
% a = 2; b = 1.5;
% xpts = a + b*exp(1i*linspace(0,2*pi, N)');
% ypts = xpts;
% A = sample_kernel(K, xpts, ypts);
% 
% svals = svd(A);
% 
% Kpnorm = 1;
% Cnorm = svals(1);
% 
% r = 25;
% % T(1)=a+b, T(-1)=a-b, T(Inf) = -T(0):
% T = @(z) sqrt(a^2-b^2)*(b-(sqrt(a^2-b^2)-a)*z)./(b+(sqrt(a^2-b^2)-a)*z); 
% R = 2*(a/b)^2 + 2*(a/b)*sqrt((a/b)^2 - 1) - 1;
% Z = R.^(-(0:r));
% zolinterperr = zeros(r+1,1);
% for n = 0:r
%     equipts = linspace(0,2*pi,n+2)';
%     equipts = equipts(1:n+1);
%     pts = T(exp(1i*equipts));
%     pls = T(R*exp(1i*equipts));
%     wts = firstBaryWeights(pts,pls);
%     zolinterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
% end
% zolinterperr = zolinterperr/svals(1);
% 
% WebbEst = Z.*Kpnorm*Cnorm/svals(1);
% 
% if discreteZn
%     aaaZ = zeros(r+1,1);
%     aaainterperr = zeros(r+1,1);
%     E = ypts;
%     F = -xpts;
%     for n = 13:20
%         [pts, pls, sigma] = aaaZolotarev(E,F,n);
%         aaaZ(n+1) = sigma;
%         wts = firstBaryWeights(pts,pls);
%         aaainterperr(n+1) = interp_error(K,xpts,ypts,pts,pls,wts);
%     end
%     aaaWebbEst = aaaZ.*Kpnorm*Cnorm/svals(1);
%     aaainterperr = aaainterperr/svals(1);
% end
% 
% svals = svals/svals(1);
% 
% %%
% 
% semilogy(0:r, WebbEst, 'color', colours(3,:))
% hold on
% semilogy(0:r, zolinterperr, 'o', 'color', colours(3,:), 'MarkerSize',3)
% if discreteZn
%     semilogy(0:r, aaaWebbEst, 'color', colours(4,:))
%     semilogy(0:r, aaainterperr, '^', 'color', colours(4,:))
% end
% semilogy(0:r, svals(1:r+1), '-*', 'color', colours(5,:))
% hold off;
% 
% xlabel('$n$ (rank)', 'Interpreter', 'latex')
% ylabel('$\|A - A_n \|_2 / \|A\|_2$', 'Interpreter', 'latex')
% axis([0,r,1e-15, 1e1])
% 
% if discreteZn
%     legend('Theorem \ref{thm:main}', 'Zolotarev interpolant', "Theorem \ref{thm:main} (discrete)", 'Discrete Zolotarev interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
% else
%     legend('Theorem \ref{thm:main}', 'Asymptotically optimal interpolant', 'Best', 'interpreter', 'latex', 'Location', 'eastoutside');
% end
% 
% set(findall(gcf,'Type','line'),'LineWidth',0.8);