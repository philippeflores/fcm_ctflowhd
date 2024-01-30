function [y,lambda] = plot_NBMblack(y,lambda,t,indVar,strLabel,varargin)

strScreen = 'full';
strLanguage = 'EN';
xlimMan = [-0.5, 4.5];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'STRLANGUAGE'
                strLanguage = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

M = size(y,2);
R = size(y{1},2);

perm = randperm(R);
y = permfacto(y,perm);
lambda = lambda(perm);

figure, pause(0.1), bigScreen(strScreen)

A = y;
for n = 1:M, for r = 1:R, A{n}(:,r) = A{n}(1:size(A{n},1),r)/max(A{n}(:,r)); end, A{n} = repmat(1-A{n}, 1, 1, 3); end

colMax = 14;
nRow = ceil(M/colMax);
nCol = ceil(M/nRow);
hRow = 4;
hLambda = 2;

subplot((hRow+1)*nRow+hLambda,nCol,(nRow*(hRow+1))*nCol+(1:(hLambda*nCol)))
stem(1:R,100*lambda,'filled','LineWidth',3,'Color','k','MarkerSize',12)
ax = gca;
ax.FontSize = 20;
if strcmpi(strLanguage,'EN')
	strLabelX = sprintf("Contribution of the %d components (in %%)",R);
elseif strcmpi(strLanguage,'FR')
	strLabelX = sprintf("Contribution des %d composantes (en %%)",R);
end
xlabel(strLabelX,'FontSize',30)
xlim([0 R+1])
xticks([1 5:5:(R-1) R])
xticklabels([1 5:5:(R-1) R])

ax1 = cell(1,M);
for row = 1:nRow
	for col = 1:nCol
		m = (row-1)*nCol+col;
		if m<M+1
			ax1{m} = subplot((hRow+1)*nRow+hLambda,nCol,(row-1)*nCol*(hRow+row-1)+(col-1)+(0:(hRow-1))*nCol+1);
			imagesc(t{indVar(m)},1:R,permute(A{m},[2 1 3]));
			colormap(1-gray)
			ax = gca; ax.FontSize = 20;

            if t{indVar(m)}(end)>100
                ecart = t{indVar(m)}(2)-t{indVar(m)}(1);
                xlim([t{indVar(m)}(1)-ecart t{indVar(m)}(end)+ecart])
            else
                xlim(xlimMan)
            end

            ax = gca; ax.FontSize = 12;
			ytickangle(0);
			if strcmpi(strLanguage,'EN')
				if m==1, ylabel("Factor Matrices",'FontSize',40), end
			elseif strcmpi(strLanguage,'FR')
				if m==1, ylabel("Matrice Facteurs",'FontSize',40), end
			end
			title(strLabel{m},'FontSize',25)
		end
	end
end

linkaxes([ax1{1} ax1{2:end}],'y')

end
