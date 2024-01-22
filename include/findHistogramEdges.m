function [a,b] = findHistogramEdges(X,varargin)

a0 = -5;
b0 = 15;
nBins = 400;
tolHist = 10^-3;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'NBINS'
                nBins = varargin{i+1};
            case 'A0'
                a0 = varargin{i+1};
            case 'B0'
                b0 = varargin{i+1};
            case 'TOLHIST'
                tolHist = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

t = linspace(a0,b0,nBins+1);
h = histcounts(X,t,'Normalization','probability');

i = 2;
j = nBins-1;

flag = 0;
while flag == 0
	if sum(h(2:i))>tolHist
		flag = 1;
	else
		i = i+1;
	end
end

flag = 0;
while flag == 0
	if sum(h(j:end-1))>tolHist
		flag = 1;
	else
		j = j-1;
	end
end

a = t(max(1,i-1));
b = t(min(nBins+1,j+1));

end