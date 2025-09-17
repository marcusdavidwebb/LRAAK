function [roots, poles, sigma] = aaaZolotarev(E,F,n)
fEF = [-ones(length(E),1); ones(length(F),1)];
[q,~,~,~,zj,fj,weights] = aaa(fEF,[E;F],'degree',n,'sign',1,'lawson',300,'damping',0.9);
tau = norm(fEF-q([E;F]),inf);
sigma = (tau/(1+sqrt(max(1-tau^2,0))))^2;
p = (1-sigma)/(1+sigma);
[~,~,roots] = prz(zj,fj+p,weights);
[~,~,poles] = prz(zj,fj-p,weights);
end