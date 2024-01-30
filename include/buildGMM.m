function [g,y,lambda] = buildGMM(X,indVar,R,t)
% This function builds a Gaussian Mixture Model from an observation matrix
% X.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

g = fitgmdist(X(:,indVar),R);

I = size(t{indVar(1)},2);
M = size(indVar,2);

y = cell(1,M);
for m = 1:M
    y{m} = zeros(I,R);
    for r = 1:R
        tempGm = gmdistribution(g.mu(r,m),g.Sigma(m,m,r));
        y{m}(:,r) = pdf(tempGm,t{indVar(m)}');
        y{m}(:,r) = y{m}(:,r)/sum(y{m}(:,r));
    end
end
lambda = g.ComponentProportion';

end

