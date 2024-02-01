function plot_Kmeans(X,indVar,t,idx,C,strLabel,varargin)

strScreen = 'full';
stepCloud = 50;
colMax = 4;
boolSeparation = 1;
xlimMan = [-0.5 4.5];
indMarg = 'all';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'STEPCLOUD'
                stepCloud = varargin{i+1};
            case 'COLMAX'
                colMax = varargin{i+1};
            case 'BOOLSEPARATION'
                boolSeparation = varargin{i+1};
            case 'INDMARG'
                indMarg = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

M = size(indVar,2);
I = size(t{indVar(1)},2);
K = max(idx);
sizeClust = zeros(K,1);
for k = 1:K
    sizeClust(k) = sum(idx==k);
end
[~,orderPlot] = sort(sizeClust,'descend');

dark2 = [55, 126, 184; 228, 26, 28; 77, 175, 74; 152, 78, 163; 255, 127, 0; 166, 86, 40]/255;
colorClust = cell(1,K);
for k = 1:K
    ind = mod(k,size(dark2,1));
    if ind==0, ind = size(dark2,1); end
    colorClust{k} = dark2(ind,:);
end

figure, bigScreen(strScreen)

if strcmpi(indMarg,'all')
    axe = cell(M);

    for m = 1:M
        axe{m,m} = subplot(M,M,(m-1)*M+m);
        ecart = t{indVar(m)}(2)-t{indVar(m)}(1);
        tempHist = histcounts(X(:,indVar(m)),'BinLimits',[t{indVar(m)}(1)-ecart/2 t{indVar(m)}(end)+ecart/2],'numBins',I,'normalization','probability');
        tempHist = tempHist/sum(tempHist);
        plot(t{indVar(m)}-ecart/2,tempHist,'k','LineWidth',3)
        ax = gca; ax.FontSize = 18;
        if m==1
            ylabel(strLabel{m},'FontSize',30,'FontWeight','bold');
            title(strLabel{m},'FontSize',30,'FontWeight','bold');
        end
        if m==M
            xlabel(strLabel{m},'FontSize',30,'FontWeight','bold');
        end

        if t{indVar(m)}(end)>100
            xlim([t{indVar(m)}(1)-ecart t{indVar(m)}(end)+ecart])
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
                    temp = (idx==k).*((1:size(X,1))'); temp(temp==0) = [];
                    scatter(X(temp(1:stepCloud:end),indVar(col)),X(temp(1:stepCloud:end),indVar(row)),10,colorClust{k},'Marker','.')
                    hold on
                end
                for k2 = 1:K
                    k = orderPlot(k2);
                    scatter(C(k,col),C(k,row),300,colorClust{k},'filled','Marker','square','LineWidth',4,'MarkerEdgeColor','k')
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
        ecart = t{indVar(indVisu1D(ind1d))}(2)-t{indVar(indVisu1D(ind1d))}(1);
        tempHist = histcounts(X(:,indVar(indVisu1D(ind1d))),'BinLimits',[t{indVar(indVisu1D(ind1d))}(1)-ecart/2 t{indVar(indVisu1D(ind1d))}(end)+ecart/2],'numBins',I,'normalization','probability');
        tempHist = tempHist/sum(tempHist);
        plot(t{indVar(indVisu1D(ind1d))}-ecart/2,tempHist,'k','LineWidth',4)
        xlabel(strLabel{indVisu1D(ind1d)},'FontSize',30,'FontWeight','bold');

        if t{indVar(indVisu1D(ind1d))}(end)>100
            xlim([t{indVisu1D(ind1d)}(1)-ecart t{indVisu1D(ind1d)}(end)+ecart])
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
            temp = (idx==k).*((1:size(X,1))'); temp(temp==0) = [];
            scatter(X(temp(1:stepCloud:end),indVar(j)),X(temp(1:stepCloud:end),indVar(i)),10,colorClust{k},'Marker','.')
            hold on
        end
        for k2 = 1:K
            k = orderPlot(k2);
            scatter(C(k,j),C(k,i),300,colorClust{k},'filled','Marker','square','LineWidth',4,'MarkerEdgeColor','k')
        end
        if t{indVar(i)}(end)>100
            ecart = t{indVar(i)}(2)-t{indVar(i)}(1);
            ylim([t{indVar(i)}(1)-ecart t{indVar(i)}(end)+ecart])
        else
            ylim(xlimMan)
        end
        if t{indVar(j)}(end)>100
            ecart = t{indVar(j)}(2)-t{indVar(j)}(1);
            xlim([t{indVar(j)}(1)-ecart t{indVar(j)}(end)+ecart])
        else
            xlim(xlimMan)
        end
    end

end
strMsgbox = strcat("{\fontsize{20}Cluster information :", sprintf("\n\nNumber of clusters : %d\n}",K));

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';

for k2 = 1:K
    k = orderPlot(k2);
    strColor = char(sprintf("%1.2f,%1.2f,%1.2f",colorClust{k}(1),colorClust{k}(2),colorClust{k}(3)));
    strK = char(strcat("{\fontsize{15}",sprintf("\nCluster #%d : \t\t",k2),"{\bf\color[rgb]{",strColor,"}",sprintf("Color"),"}",sprintf("\nProportion : "),"{\bf",sprintf("%2.1f%%",100*sizeClust(k)/size(X,1)),"}",sprintf("\n"),"}"));
    strMsgbox = strcat(strMsgbox, strK);
end

h = msgbox(strMsgbox,"Cluster Information",CreateStruct);
temp = h.OuterPosition;

set(h,"OuterPosition",[0 4000 temp(3) temp(4)])

end

