function dataMarg = computeMarginals(X,indVar,t,I,calT)

M = size(indVar,2);

dataMarg.M = M;
dataMarg.I = I;
T = size(calT,2);
dataMarg.T = T;
dataMarg.calT = calT;
dataMarg.indVar = indVar;

marg = cell(1,T);

if T>0, fprintf("Computing the marginals..."), end
for j = 1:T
	if j>1, fprintf(repmat('\b',1,length(strHisto))), end
	strHisto = sprintf(' (%d over %d)',j,T);
	fprintf(strHisto)
	marg{j} = histnd(X(:,calT{j}(1)),X(:,calT{j}(2)),X(:,calT{j}(3)),t{calT{j}(1)},t{calT{j}(2)},t{calT{j}(3)});
	marg{j} = marg{j}/sum(sum(sum(marg{j})));
end
if T>0, fprintf("\tDone.\n"), end

dataMarg.marg = marg;

end
