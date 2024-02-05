function plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar,varargin)
% plot_Marg : plots CTFlowHD results with a marginal visualization. For
% more information on this matter, please refer to the thesis manuscript of
% the author.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

strScreen = 'full';
stepCloud = 50;
colMax = 4;
boolSeparation = 1;
xlimMan = [-0.5 4.5];
indGroup = 1:size(compGroup,2);
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
            case 'INDGROUP'
                indGroup = varargin{i+1};
            case 'INDMARG'
                indMarg = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

compGroupLoc = compGroup(indGroup);

M = size(y,2);
I = size(y{1},1);

temp = zeros(size(compGroupLoc,2),1);
for ng = 1:size(temp,1)
    temp(ng) = sum(lambda(compGroupLoc{ng}.indGroup));
end
[~,indNg] = sort(temp,'descend');

figure, bigScreen(strScreen)

if strcmpi(indMarg,'all')
    axe = cell(M);

    for m = 1:M
        axe{m,m} = subplot(M,M,(m-1)*M+m);
        ecart = t{indVar(m)}(2)-t{indVar(m)}(1);
        tempHist = histcounts(X(:,indVar(m)),'BinLimits',[t{indVar(m)}(1)-ecart/2 t{indVar(m)}(end)+ecart/2],'numBins',I,'normalization','probability');
        tempHist = tempHist/sum(tempHist);

        for ng2 = 1:size(compGroupLoc,2)
            ng = indNg(ng2);
            temp = zeros(I,1);
            for r = 1:compGroupLoc{ng}.sizeGroup
                temp = temp+lambda(compGroupLoc{ng}.indGroup(r))*y{m}(:,compGroupLoc{ng}.indGroup(r));
            end
            plot(t{indVar(m)},temp,'LineWidth',4,'Color',compGroupLoc{ng}.colorGroup)
            hold on
        end
        plot(t{indVar(m)}-ecart/2,tempHist,'k--','LineWidth',2)
        if m==1
            ylabel(strLabel{m},'FontSize',35,'FontWeight','bold');
            title(strLabel{m},'FontSize',35,'FontWeight','bold');
        end
        if m==M
            xlabel(strLabel{m},'FontSize',35,'FontWeight','bold');
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
                scatterColor(X(:,indVar(col)),X(:,indVar(row)),'stepCloud',stepCloud,'newFig','off')
                if t{indVar(row)}(end)>100
                    ecart = t{indVar(row)}(2)-t{indVar(row)}(1);
                    ylim([t{indVar(row)}(1)-ecart t{indVar(row)}(end)+ecart])
                else
                    ylim(xlimMan)
                end
                if t{indVar(col)}(end)>100
                    ecart = t{indVar(col)}(2)-t{indVar(col)}(1);
                    xlim([t{indVar(col)}(1)-ecart t{indVar(col)}(end)+ecart])
                else
                    xlim(xlimMan)
                end
                hold on
                for ng2 = 1:size(compGroupLoc,2)
                    ng = indNg(ng2);
                    temp = zeros(I);
                    for r = 1:compGroupLoc{ng}.sizeGroup
                        temp = temp+lambda(compGroupLoc{ng}.indGroup(r))*y{row}(:,compGroupLoc{ng}.indGroup(r))*y{col}(:,compGroupLoc{ng}.indGroup(r))';
                    end
                    temp = padarray(temp,[1 1],0,'both');
                    [~,cLine] = contour(axe{row,col},[2*t{indVar(col)}(1)-t{indVar(col)}(2) t{indVar(col)} 2*t{indVar(col)}(end)-t{indVar(col)}(end-1)],[2*t{indVar(row)}(1)-t{indVar(row)}(2) t{indVar(row)} 2*t{indVar(row)}(end)-t{indVar(row)}(end-1)],temp*I,linspace(0,sum(sum(temp))/(10+90*sum(lambda(compGroup{ng}.indGroup))),2));
                    cLine.LineColor = compGroupLoc{ng}.colorGroup;
                    cLine.LineWidth = 4;
                    ax = gca; ax.FontSize = 18;
                    if row==M
                        xlabel(strLabel{col},'FontSize',35,'FontWeight','bold')
                    end
                    if col==1
                        ylabel(strLabel{row},'FontSize',35,'FontWeight','bold')
                    end
                    if row==1
                        title(strLabel{col},'FontSize',35,'FontWeight','bold')
                    end
                end
            end
        end
    end

    for m = 1:M
        for row = 1:M
            if row~=m
                linkaxes([axe{m,m} axe{row,m}],'x')
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

        for ng2 = 1:size(compGroupLoc,2)
            ng = indNg(ng2);
            temp = zeros(I,1);
            for r = 1:compGroupLoc{ng}.sizeGroup
                temp = temp+lambda(compGroupLoc{ng}.indGroup(r))*y{indVisu1D(ind1d)}(:,compGroupLoc{ng}.indGroup(r));
            end
            plot(t{indVar(indVisu1D(ind1d))},temp,'LineWidth',4,'Color',compGroupLoc{ng}.colorGroup)
            hold on
        end
        plot(t{indVar(indVisu1D(ind1d))}-ecart/2,tempHist,'k--','LineWidth',2)
        xlabel(strLabel{indVisu1D(ind1d)},'FontSize',35,'FontWeight','bold');

        if t{indVar(indVisu1D(ind1d))}(end)>100
            xlim([t{indVisu1D(ind1d)}(1)-ecart t{indVisu1D(ind1d)}(end)+ecart])
        else
            xlim(xlimMan)
        end

    end

    ax1 = cell(size(indVisu2D,1));
    for ind2d = 1:size(indVisu2D,1)
        j = indVisu2D(ind2d,1);
        i = indVisu2D(ind2d,2);
        ax1{ind2d} = subplot(nRow,nCol,indSubPlot(ind2d+size(indVisu1D,2)));
        scatterColor(X(:,indVar(j)),X(:,indVar(i)),'stepCloud',stepCloud,'newFig','off')
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
        hold on
    end

    for ind2d = 1:size(indVisu2D,1)
        subplot(nRow,nCol,indSubPlot(ind2d+size(indVisu1D,2)))
        for ng2 = 1:size(compGroupLoc,2)
            ng = indNg(ng2);
            pause(0.1)
            j = indVisu2D(ind2d,1);
            i = indVisu2D(ind2d,2);
            temp = zeros(I);
            for r = 1:compGroupLoc{ng}.sizeGroup
                temp = temp+lambda(compGroupLoc{ng}.indGroup(r))*y{i}(:,compGroupLoc{ng}.indGroup(r))*y{j}(:,compGroupLoc{ng}.indGroup(r))';
            end
            temp = padarray(temp,[1 1],0,'both');
            [~,cLine] = contour(ax1{ind2d},[2*t{indVar(j)}(1)-t{indVar(j)}(2) t{indVar(j)} 2*t{indVar(j)}(end)-t{indVar(j)}(end-1)],[2*t{indVar(i)}(1)-t{indVar(i)}(2) t{indVar(i)} 2*t{indVar(i)}(end)-t{indVar(i)}(end-1)],temp*I,linspace(0,sum(sum(temp))/(10+90*sum(lambda(compGroup{ng}.indGroup))),2));
            cLine.LineColor = compGroupLoc{ng}.colorGroup;
            cLine.LineWidth = 4;
            ax = gca; ax.FontSize = 18;
            xlabel(strLabel{j},'FontSize',35,'FontWeight','bold');
            ylabel(strLabel{i},'FontSize',35,'FontWeight','bold');
        end
    end
end

end
