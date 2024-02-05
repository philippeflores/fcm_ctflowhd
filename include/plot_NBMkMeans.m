function plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,varargin)
% plot_NBMkMeans : plots CTFlowHD results with a K-means visualization. For
% more information on this matter, please refer to the thesis manuscript of
% the author.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

strScreen = 'full';
colMax = 4;
boolSeparation = 1;
xlimMan = [-0.5 4.5];
indMarg = 'all';
strCentroid = 'esp';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'COLMAX'
                colMax = varargin{i+1};
            case 'BOOLSEPARATION'
                boolSeparation = varargin{i+1};
            case 'STRCENTROID'
                strCentroid = varargin{i+1};
            case 'INDMARG'
                indMarg = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

M = size(y,2);
I = size(y{1},1);
R = size(y{1},2);

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

[idx,C] = kmeans(matFeature,K);

sizeClust = zeros(K,1);
propClust = zeros(K,1);
for k = 1:K
    sizeClust(k) = sum(idx==k);
    propClust(k) = sum(lambda(idx==k));
end
[~,orderPlot] = sort(sizeClust,'descend');

dark2 = [55, 126, 184; 228, 26, 28; 77, 175, 74; 152, 78, 163; 255, 127, 0; 166, 86, 40]/255; nColor = size(dark2,1);
colorClust = cell(1,K);
for k2 = 1:K
    k = orderPlot(k2);
    ind = mod(k,size(dark2,1));
    if ind==0, ind = size(dark2,1); end
    colorClust{k} = dark2(ind,:);
end

shapes = {'o','v','square','^','<','diamond','>'}; nShapes = size(shapes,2);
shapeClust = cell(1,K);
for k2 = 1:K
    k = orderPlot(k2);
	if k<=nShapes*nColor
		shapeClust{k} = shapes{ceil(k/nColor)};
	else
		shapeClust{k} = shapes{1};
	end
end

figure, bigScreen(strScreen)

if strcmpi(indMarg,'all')
    axe = cell(M);

    for m = 1:M
        axe{m,m} = subplot(M,M,(m-1)*M+m);

        temp = zeros(I,1);
        for r = 1:R
            temp = temp+lambda(r)*y{m}(:,r);
        end
        temp = temp/sum(temp);

        plot(t{indVar(m)},temp,'k','LineWidth',3)
        
        hold on
        for k2 = 1:K
            k = orderPlot(k2);
            temp = zeros(I,1);
            for r = 1:R
                if idx(r)==k
                    temp = temp+lambda(r)*y{m}(:,r);
                end
            end
            plot(t{indVar(m)},temp,'--','Color',colorClust{k},'LineWidth',2)
        end

        ax = gca; ax.FontSize = 18;
        if m==1
            ylabel(strLabel{m},'FontSize',30,'FontWeight','bold');
            title(strLabel{m},'FontSize',30,'FontWeight','bold');
        end
        if m==M
            xlabel(strLabel{m},'FontSize',30,'FontWeight','bold');
        end

        if t{indVar(m)}(end)>100
            xlim([t{indVar(m)}(1) t{indVar(m)}(end)])
        else
            xlim(xlimMan)
        end
    end

    for row = 1:M
        for col = 1:M
            if row~=col
                axe{row,col} = subplot(M,M,(row-1)*M+col);
                for k2 = 1:K
                    k = orderPlot(k2);
                    temp = (idx==k).*((1:R)'); temp(temp==0) = [];
                    scatter(matFeature(temp,col),matFeature(temp,row),250*R*lambda(temp),colorClust{k},'filled','Marker',shapeClust{k},'MarkerFaceAlpha',0.2,'MarkerEdgeColor','k')
                    hold on
                end
                for k2 = 1:K
                    k = orderPlot(k2);
                    scatter(C(k,col),C(k,row),15*R,colorClust{k},'filled','Marker',shapeClust{k},'LineWidth',2,'MarkerEdgeColor','k')
                end

                if t{indVar(row)}(end)>100
                    ylim([t{indVar(row)}(1) t{indVar(row)}(end)])
                else
                    ylim(xlimMan)
                end

                ax = gca; ax.FontSize = 18;
                if row==M
                    xlabel(strLabel{col},'FontSize',30,'FontWeight','bold')
                end
                if col==1
                    ylabel(strLabel{row},'FontSize',30,'FontWeight','bold')
                end
                if row==1
                    title(strLabel{col},'FontSize',30,'FontWeight','bold')
                end
                linkaxes([axe{col,col} axe{row,col}],'x')
            end
        end
    end

else

    indVisu1D = zeros(1,size(indMarg,2));
    indVisu2D = zeros(size(indMarg,2),2);
    for i = 1:size(indMarg,2)
        if size(indMarg{i},2)==1
            indVisu1D(i) = indMarg{i};
        elseif size(indMarg{i},2)==2
            indVisu2D(i,:) = indMarg{i};
        end
    end
    indVisu2D(indVisu1D>0,:) = [];
    indVisu1D(indVisu1D==0) = [];

    nPlot = size(indVisu1D,2)+size(indVisu2D,1);

    if boolSeparation==0
        nRow = ceil(nPlot/colMax);
        nCol = min([colMax nPlot]);
        indSubPlot = 1:nPlot;
    else
        if size(indVisu1D,2)>size(indVisu2D,1)
            nCol = size(indVisu1D,2);
            nRow = 2;
            indSubPlot = 1:nPlot;
        else
            nCol = min([size(indVisu2D,1) colMax]);
            nRow = ceil(size(indVisu1D,2)/nCol)+ceil(size(indVisu2D,1)/nCol);
            indSubPlot = [1:size(indVisu1D,2) ceil(size(indVisu1D,2)/nCol)*nCol+(1:size(indVisu2D,1))];
        end
    end

    for ind1d = 1:size(indVisu1D,2)
        subplot(nRow,nCol,indSubPlot(ind1d))
        
        temp = zeros(I,1);
        for r = 1:R
            temp = temp+lambda(r)*y{ind1d}(:,r);
        end
        temp = temp/sum(temp);
        plot(t{indVar(ind1d)},temp,'k','LineWidth',3)
        hold on
        for k2 = 1:K
            k = orderPlot(k2);
            temp = zeros(I,1);
            for r = 1:R
                if idx(r)==k
                    temp = temp+lambda(r)*y{ind1d}(:,r);
                end
            end
            plot(t{indVar(ind1d)},temp,'--','Color',colorClust{k},'LineWidth',2)
        end

        xlabel(strLabel{indVisu1D(ind1d)},'FontSize',30,'FontWeight','bold');

        if t{indVar(indVisu1D(ind1d))}(end)>100
            xlim([t{indVisu1D(ind1d)}(1) t{indVisu1D(ind1d)}(end)])
        else
            xlim(xlimMan)
        end

    end


    for ind2d = 1:size(indVisu2D,1)
        j = indVisu2D(ind2d,1);
        i = indVisu2D(ind2d,2);
        subplot(nRow,nCol,indSubPlot(ind2d+size(indVisu1D,2)))
        for k2 = 1:K
            k = orderPlot(k2);
            temp = (idx==k).*((1:R)'); temp(temp==0) = [];
            scatter(matFeature(temp,j),matFeature(temp,i),250*R*lambda(temp),colorClust{k},'filled','Marker',shapeClust{k},'MarkerFaceAlpha',0.2,'MarkerEdgeColor','k')
            hold on
        end
        for k2 = 1:K
            k = orderPlot(k2);
            scatter(C(k,j),C(k,i),15*R,colorClust{k},'filled','Marker',shapeClust{k},'LineWidth',2,'MarkerEdgeColor','k')
        end
        if t{indVar(i)}(end)>100
            ylim([t{indVar(i)}(1) t{indVar(i)}(end)])
        else
            ylim(xlimMan)
        end
        if t{indVar(j)}(end)>100
            xlim([t{indVar(j)}(1) t{indVar(j)}(end)])
        else
            xlim(xlimMan)
        end
        xlabel(strLabel{j},'FontSize',35,'FontWeight','bold');
        ylabel(strLabel{i},'FontSize',35,'FontWeight','bold');
    end

end

strMsgbox = strcat("{\fontsize{20}Cluster information :", sprintf("\n\nNumber of clusters : %d\n}",K));

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';

for k2 = 1:K
    k = orderPlot(k2);
    strColor = char(sprintf("%1.2f,%1.2f,%1.2f",colorClust{k}(1),colorClust{k}(2),colorClust{k}(3)));
    strK = char(strcat("{\fontsize{15}",sprintf("\nCluster #%d : \t\t",k2),"{\bf\color[rgb]{",strColor,"}",sprintf("Color\t\t%s",shapeClust{k}),"}",sprintf("\nProportion : "),"{\bf",sprintf("%2.1f%%",100*propClust(k)),"}",sprintf("\n"),"}"));
    strMsgbox = strcat(strMsgbox, strK);
end

h = msgbox(strMsgbox,"Cluster Information",CreateStruct);
temp = h.OuterPosition;

set(h,"OuterPosition",[0 4000 temp(3) temp(4)])

end

