function Y = plot_NBMtSNE(y,lambda,compGroup,t,strLabel,indVar,varargin)

colMax = 4;
perplexity = -1;
strScreen = 'halfR';
strCentroid = 'esp';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'COLMAX'
                colMax = varargin{i+1};
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'STRCENTROID'
                strCentroid = varargin{i+1};
            case 'PERPLEXITY'
                perplexity = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

R = size(y{1},2);
if perplexity<0, perplexity = max([2 R/4]); end
nGroup = size(compGroup,2);
M = size(y,2);
I = size(y{1},1);


k = 2;
colMax = colMax+2;
if M>10
    nRow = ceil((M-(colMax*2-4))/colMax)+2;
else
    nRow = ceil((M+k*k)/colMax);
end
nCol = min(colMax,k+ceil(M/nRow));
indK = zeros(k*k,1); for ind = 1:k*k, modind = mod(ind,k); modind(modind==0)=k; indK(ind) = (floor((ind-1)/k)*nCol)+modind; end
indRem = 1:(nRow*nCol); indRem(indK)=[];

matFeature = zeros(R,M);
if strcmpi(strCentroid,'esp')
    for r = 1:R
        for m = 1:M
            matFeature(r,m) = t{indVar(m)}*y{m}(:,r);
        end
    end
elseif strcmpi(strCentroid,'max')
	for r = 1:R
		for m = 1:M
			[~,matFeature(r,m)] = max(y{m}(:,r));
            matFeature(r,m) = t{indVar(m)}(matFeature(r,m));
		end
	end
end

figure
bigScreen(strScreen)

Y = tsne(matFeature,"Perplexity",perplexity);

ax = subplot(nRow,nCol,indK);
hold on
for ng = 1:nGroup
    structGroup = compGroup{ng};
    [~,orderPlot] = sort(lambda(structGroup.indGroup),'descend');
    for r = 1:size(structGroup.indGroup,2)
        a = plot(Y(structGroup.indGroup(orderPlot(r)),1),Y(structGroup.indGroup(orderPlot(r)),2),'Marker',compGroup{ng}.shapeGroup,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',structGroup.colorGroup,'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*50*R*k*k/(nRow*nCol));
        row11 = dataTipTextRow('R',structGroup.indGroup(orderPlot(r)));
        a.DataTipTemplate.DataTipRows(1) = row11;
        row12 = dataTipTextRow('R_{first}',structGroup.indGroup(1));
        a.DataTipTemplate.DataTipRows(2) = row12;
        row13 = dataTipTextRow('R_{last}',structGroup.indGroup(end));
        a.DataTipTemplate.DataTipRows(end+1) = row13;
        row2 = dataTipTextRow('\lambda (%)',100*lambda(structGroup.indGroup(orderPlot(r))));
        a.DataTipTemplate.DataTipRows(end+1) = row2;
        row3 = dataTipTextRow('\lambda_{Total} (%)',100*sum(lambda(structGroup.indGroup)));
        a.DataTipTemplate.DataTipRows(end+1) = row3;
    end
    ax.FontSize = 15;
    title("Hierarchical clustering",'FontSize',35)
end

axExp = cell(1,M);
cmapJet = colormap([50,136,189; 102,194,165; 171,221,164; 230,245,152; 255,255,191; 254,224,139; 253,174,97; 244,109,67; 213,62,79]/255);
for m = 1:M
    posEsp = round(sum(y{m}.*(repmat(1:I,R,1)')));
    for ng = 1:nGroup
        structGroup = compGroup{ng};
        [~,orderPlot] = sort(lambda(structGroup.indGroup),'descend');
        axExp{m} = subplot(nRow,nCol,indRem(m));
        hold on
        for r = 1:size(structGroup.indGroup,2)
            a = plot(Y(structGroup.indGroup(orderPlot(r)),1),Y(structGroup.indGroup(orderPlot(r)),2),'Marker',compGroup{ng}.shapeGroup,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cmapJet(ceil(posEsp(structGroup.indGroup(orderPlot(r)))*size(cmapJet,1)/I),:),'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*50*R*k/(nRow*nCol));
            row11 = dataTipTextRow('R',structGroup.indGroup(orderPlot(r)));
            a.DataTipTemplate.DataTipRows(1) = row11;
            row12 = dataTipTextRow('R_{first}',structGroup.indGroup(1));
            a.DataTipTemplate.DataTipRows(2) = row12;
            row13 = dataTipTextRow('R_{last}',structGroup.indGroup(end));
            a.DataTipTemplate.DataTipRows(end+1) = row13;
            row2 = dataTipTextRow('\lambda (%)',100*lambda(structGroup.indGroup(orderPlot(r))));
            a.DataTipTemplate.DataTipRows(end+1) = row2;
            row3 = dataTipTextRow('\lambda_{Total} (%)',100*sum(lambda(structGroup.indGroup)));
            a.DataTipTemplate.DataTipRows(end+1) = row3;
        end
    end
    if t{indVar(m)}(end)>100
        clim([t{indVar(m)}(1) t{indVar(m)}(end)])
    else
        clim([-0.5 4.5])
    end
    cb = colorbar(); cb.Label.Position(1) = 4;
    if t{indVar(m)}(end)<100
        ylabel(cb,'Fluorescence','FontSize',18,'Rotation',270)
    end
    axExp{m}.FontSize = 15;
    title(sprintf("%s",strLabel{m}),'FontSize',35)

end
linkaxes([ax axExp{:}],'xy')


end