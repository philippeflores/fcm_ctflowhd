function plot_GMM(g,X,indVar,t,strLabel,varargin)

stepCloud = 10;
indPlot = {};
colMax = -1;
strScreen = 'full';
strXlim = 'support';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STEPCLOUD'
                stepCloud = varargin{i+1};
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'INDPLOT'
                indPlot = varargin{i+1};
            case 'COLMAX'
                colMax = varargin{i+1};
            case 'STRXLIM'
                strXlim = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

R = g.NumComponents;
M = g.NumVariables;
I = size(t{indVar(1)},2);
if colMax==-1, colMax=M; end

colorR = cell(1,R);

figure
bigScreen(strScreen)

if isempty(indPlot)==1

    ax2D = cell(M,M);
    ax1D = cell(1,M);

    for m1 = 1:size(indVar,2)
        for m2 = 1:size(indVar,2)
            if m1~=m2
                ax2D{m2,m1} = subplot(M,M,(m1-1)*M+m2);
                plot(X(1:stepCloud:end,indVar(m2)),X(1:stepCloud:end,indVar(m1)),'o', ...
                    'LineStyle','none','Color',0.4*[1 1 1],'MarkerSize',3)
                if strcmpi(strXlim,'auto')==0
                    if strcmpi(strXlim,'support')
                        xlim([t{indVar(m2)}(1) t{indVar(m2)}(end)])
                        ylim([t{indVar(m1)}(1) t{indVar(m1)}(end)])
                    elseif strcmpi(strXlim,'fluo')
                        xlim([-0.5 4.5]), ylim([-0.5 4.5])
                    end
                end
                hold on
                for r = 1:R
                    tempGm = gmdistribution(g.mu(r,[m2 m1]),g.Sigma([m2 m1],[m2 m1],r));
                    tempgmPDF = @(x,y) arrayfun(@(x0,y0) pdf(tempGm,[x0 y0]),x,y);
                    if strcmpi(strXlim,'fluo')
                        tempSupport = [-0.5 4.5];
                    else
                        tempSupport = [t{indVar(m2)}(1) t{indVar(m2)}(end) t{indVar(m1)}(1) t{indVar(m1)}(end)];
                    end
                    fcontour(tempgmPDF,tempSupport,'LineColor','k')
                    tempScat = scatter(g.mu(r,m2),g.mu(r,m1),150,'filled');
                    colorR{r} = tempScat.CData;
                end
                if m1==1
                    title(strLabel{m2},'FontSize',25)
                end
                if m2==1
                    ylabel(strLabel{m1},'FontWeight','bold','FontSize',25)
                end
            end
        end
    end
    for m = 1:size(indVar,2)
        ax1D{m} = subplot(M,M,(m-1)*M+m);
        htemp = histcounts(X(:,indVar(m)),I,'BinLimits',t{indVar(m)}([1 I]),'Normalization','pdf');
        plot(t{indVar(m)},htemp,'k-','LineWidth',2)
        if strcmpi(strXlim,'auto')==0
            if strcmpi(strXlim,'support')
                xlim([t{indVar(m)}(1) t{indVar(m)}(end)])
            elseif strcmpi(strXlim,'fluo')
                xlim([-0.5 4.5])
            end
        end
        hold on
        for r = 1:R
            tempGm = gmdistribution(g.mu(r,m),g.Sigma(m,m,r));
            tempgmPDF = @(x) arrayfun(@(x0) g.ComponentProportion(r)*pdf(tempGm,x0),x);
            if strcmpi(strXlim,'fluo')
                tempSupport = [-0.5 4.5];
            else
                tempSupport = [t{indVar(m)}(1) t{indVar(m)}(end)];
            end
            fplot(tempgmPDF,tempSupport,':','Color',colorR{r},'LineWidth',2)
        end
        if m==1
            title(strLabel{m},'FontSize',25)
            ylabel(strLabel{m},'FontWeight','bold','FontSize',25)
        end
    end

    for m = 1:M
        temp = [];
        for m2 = 1:M
            if m2~=m
                temp = [temp ax2D{m,m2}];
            end
        end
        linkaxes([ax1D{m} temp],'x')
    end


    for m = 1:M
        temp = [];
        for m1 = 1:M
            if m1~=m
                temp = [temp ax2D{m1,m}];
            end
        end
        linkaxes(temp,'y')
    end

else

    nMarg = size(indPlot,2);
    nMarg1 = 0; nMarg2 = 0; temp = ones(nMarg,1); for n = 1:nMarg, if size(indPlot{n},2)==1, nMarg1 = nMarg1+1; else, nMarg2 = nMarg2+1; temp(n) = 2; end, end
    nC = max([min([colMax nMarg1]) min([colMax nMarg2])]);
    nRow = ceil(nMarg1/nC)+ceil(nMarg2/nC);
    indMarg1 = cell2mat(indPlot(temp==1));
    indMarg2 = cell2mat(indPlot(temp==2)');
    for n = 1:nMarg1
        m = indMarg1(n);
        subplot(nRow,nC,n)
        htemp = histcounts(X(:,indVar(m)),I,'BinLimits',t{indVar(m)}([1 I]),'Normalization','pdf');
        plot(t{indVar(m)},htemp,'k-','LineWidth',2)
        if strcmpi(strXlim,'auto')==0
            if strcmpi(strXlim,'support')
                xlim([t{indVar(m)}(1) t{indVar(m)}(end)])
            elseif strcmpi(strXlim,'fluo')
                xlim([-0.5 4.5])
            end
        end
        hold on
        for r = 1:R
            tempGm = gmdistribution(g.mu(r,m),g.Sigma(m,m,r));
            tempgmPDF = @(x) arrayfun(@(x0) g.ComponentProportion(r)*pdf(tempGm,x0),x);
            if strcmpi(strXlim,'fluo')
                tempSupport = [-0.5 4.5];
            else
                tempSupport = [t{indVar(m)}(1) t{indVar(m)}(end)];
            end
            tempPlot = fplot(tempgmPDF,tempSupport,':','LineWidth',2);
            colorR{r} = tempPlot.Color;
        end
        title(strLabel{m},'FontSize',25)
    end
    offset = ceil(nMarg1/nC)*nC;
    for n = 1:nMarg2
        subplot(nRow,nC,offset+n)
        m1 = indMarg2(n,1);
        m2 = indMarg2(n,2);
        plot(X(1:stepCloud:end,indVar(m2)),X(1:stepCloud:end,indVar(m1)),'.', ...
            'LineStyle','none','Color',0.7*[1 1 1],'MarkerSize',7)
        if strcmpi(strXlim,'auto')==0
            if strcmpi(strXlim,'support')
                xlim([t{indVar(m2)}(1) t{indVar(m2)}(end)])
                ylim([t{indVar(m1)}(1) t{indVar(m1)}(end)])
            elseif strcmpi(strXlim,'fluo')
                xlim([-0.5 4.5]), ylim([-0.5 4.5])
            end
        end
        hold on
        for r = 1:R
            tempGm = gmdistribution(g.mu(r,[m2 m1]),g.Sigma([m2 m1],[m2 m1],r));
            tempgmPDF = @(x,y) arrayfun(@(x0,y0) pdf(tempGm,[x0 y0]),x,y);

            if isempty(colorR{r})==0
                scatter(g.mu(r,m2),g.mu(r,m1),250,colorR{r},'filled');
                if strcmpi(strXlim,'fluo')
                    tempSupport = [-0.5 4.5];
                else
                    tempSupport = [t{indVar(m2)}(1) t{indVar(m2)}(end) t{indVar(m1)}(1) t{indVar(m1)}(end)];
                end
                fcontour(tempgmPDF,tempSupport,'LineColor',colorR{r})
            else
                temp = scatter(g.mu(r,m2),g.mu(r,m1),250,'filled');
                colorR{r} = temp.CData;
                if strcmpi(strXlim,'fluo')
                    tempSupport = [-0.5 4.5];
                else
                    tempSupport = [t{indVar(m2)}(1) t{indVar(m2)}(end) t{indVar(m1)}(1) t{indVar(m1)}(end)];
                end
                fcontour(tempgmPDF,tempSupport,'LineColor',colorR{r},'LineWidth',2)
            end
        end
        ylabel(strLabel{m1},'FontWeight','bold','FontSize',25)
        xlabel(strLabel{m2},'FontWeight','bold','FontSize',25)
    end

end

end