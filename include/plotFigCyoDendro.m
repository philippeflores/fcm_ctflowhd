function [y,lambda,compGroup,distR,Z,H] = plotFigCytoDendro(y,lambda,t,indVar,strFluorescence,colorThresh,varargin)

strScreen = 'b';
strLanguage = 'EN';
xlimMan = [-0.5, 4.5];
boolRegroup = 0;
strDist = 'esp';
strLink = 'complete';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'STRLANGUAGE'
                strLanguage = varargin{i+1};
            case 'XLIMMAN'
                xlimMan = varargin{i+1};
            case 'BOOLREGROUP'
                boolRegroup = varargin{i+1};
            case 'STRDIST'
                strDist = varargin{i+1};
            case 'STRLINK'
                strLink = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

N = size(y,2);
R = size(y{1},2);
nHist = size(y{1},1);
boolNorm = 1;

figDendro = figure; figDendro = figDendro.Number;

pause(0.1)

if strcmpi(strScreen,'b')
	vScreen =  [-1.43, 0, 1.43,1.33];
elseif strcmpi(strScreen,'s')
	vScreen = [0, 0, 1,1];
end

B = y;
for n = 1:N, for r = 1:R, B{n}(:,r) = B{n}(1:size(B{n},1),r)/max(B{n}(:,r)); end, B{n} = repmat(1-B{n}, 1, 1, 3); end

Nmax = 14;
nRow = ceil(N/Nmax);
nCol = ceil(N/nRow);
hRow = 4;
hLambda = 2;

set(gcf, 'Units', 'Normalized', 'OuterPosition', vScreen);

if strcmpi(strDist,'max')
	distR = zeros(R,N);
	for r = 1:R
		for n = 1:N
			[~,distR(r,n)] = max(y{n}(:,r));
		end
	end
elseif strcmpi(strDist,'esp2')
	distR = zeros(R,N);
	for r = 1:R
		for n = 1:N
			distR(r,n) = (sum((1:nHist)'.*y{n}(:,r)));
		end
	end
elseif strcmpi(strDist,'esp')
	distR = zeros(R,N);
	for r = 1:R
		for n = 1:N
			distR(r,n) = (t{indVar(n)}*y{n}(:,r))/max(t{indVar(n)});
		end
	end
elseif strcmpi(strDist,'corr')
	distR = ones(R)-eye(R);
	for r = 1:R
		for s = 1:R
			if r~=s
				for n = 1:N
					if boolNorm
						tempR = y{n}(:,r)/norm(y{n}(:,r));
						tempS = y{n}(:,s)/norm(y{n}(:,s));
					else
						tempR = y{n}(:,r);
						tempS = y{n}(:,s);
					end
					distR(r,s) = distR(r,s)*(dot(tempR,tempS));
				end
				distR(r,s) = 1-(distR(r,s))^(1/(2^N));
			end
		end
	end
else
	error('Chosen distance unavailable...\n');
end

subplot(hRow*nRow+hLambda+hLambda+nRow+1,nCol,(nRow*hRow+nRow)*nCol+[1:2*nCol])
if strcmpi(strDist,'corr')
	Z = linkage(squareform(distR),strLink);
else
	Z = linkage(distR,strLink);
end

Z(:,3) = Z(:,3)/max(Z(:,3));
[H,~,outperm,palColor] = myDendrogram(Z,R,0,'ColorThreshold',colorThresh);
hold on
plot(xlim(),[1 1]*colorThresh,'--r','LineWidth',4)
xticklabels(string(1:R))
set(H,'Linewidth',5)
ax = gca; ax.FontSize = 15;
% xtickangle(70);
if strcmpi(strLanguage,'EN')
	ylabel("Dendrogram",'FontSize',30);
elseif strcmpi(strLanguage,'FR')
	ylabel("Dendrogramme",'FontSize',40);
end
ylim([0 1.05])

lambda = lambda(outperm);
y = permfacto(y,outperm);

for r = 1:R-1
	if Z(r,1)<R+1
		flag = 0; s = 1;
		while flag==0
			if Z(r,1)==outperm(s)
				flag = 1;
				Z(r,1) = s;
			else
				s = s+1;
			end
		end
	end
end

for r = 1:R-1
	if Z(r,2)<R+1
		flag = 0; s = 1;
		while flag==0
			if Z(r,2)==outperm(s)
				flag = 1;
				Z(r,2) = s;
			else
				s = s+1;
			end
		end
	end
end

compGroup = changeFigColorDendro(figDendro,palColor,y,lambda,Z,t,strFluorescence,indVar,strLanguage,H,colorThresh,xlimMan,boolRegroup);

end






function compGroup = changeFigColorDendro(figNumber,palColor,y,lambda,Z,t,strFluorescence,indVar,strLanguage,H,colorThresh,xlimMan,boolRegroup)

% On récupère les variables
R = size(y{1},2);
I = size(y{1},1);
M = size(y,2);

% On fixe les variables pour l'affichage
Nmax = 14;
nRow = ceil(M/Nmax);
nCol = ceil(M/nRow);
hRow = 4;
hLambda = 2;

% Recherche du nombre de groupes
palColorUnique = unique(palColor,"rows"); palColorUnique = palColorUnique(sum(palColorUnique,2)~=0,:);
nGroup = size(palColorUnique,1);
for r = 1:R
	indZ = mod(find(Z==r),R-1);
	if indZ == 0, indZ = R-1; end
	if length(indZ)>1, indZ = indZ(1); end
	if sum(palColor(indZ,:))==0, nGroup = nGroup+1; end
end

% Assignation des composantes à un groupe
compGroup = cell(1,nGroup);
countBlack = 1;
for ng = 1:nGroup
	compGroup{ng}.indGroup = [];
	compGroup{ng}.sizeGroup = 0;
end
for r = 1:R
	indZ = mod(find(Z==r),R-1);
	if indZ == 0, indZ = R-1; end
	if length(indZ)>1, indZ = indZ(1); end
	if sum(palColor(indZ,:)==0)==3
		compGroup{size(palColorUnique,1)+countBlack}.indGroup = r;
		compGroup{size(palColorUnique,1)+countBlack}.sizeGroup = 1;
		countBlack = countBlack+1;
	else
		ng = nonzeros((sum(palColorUnique==palColor(indZ,:),2)==3).*(1:size(palColorUnique,1))');
		compGroup{ng}.indGroup = [compGroup{ng}.indGroup r];
	end
end
for ng = 1:size(palColorUnique,1), compGroup{ng}.sizeGroup = size(compGroup{ng}.indGroup,2); end

% Tri des groupes
indPermute = zeros(nGroup,1);
for ng = 1:nGroup, indPermute(ng) = compGroup{ng}.indGroup(1); end
[~,indPermute] = sort(indPermute);
compGroup = compGroup(indPermute);

% Chargement des couleurs et dégradés
dark2 = [55, 126, 184; 228, 26, 28; 77, 175, 74; 152, 78, 163; 255, 127, 0; 166, 86, 40]/255;
nColor = size(dark2,1);
shades = cell(1,nColor);
shades{1} = [55, 126, 184; 75, 134, 184; 93, 143, 184; 110, 151, 184]/255;
shades{2} = [228, 26, 28; 233, 73, 69; 232, 106, 104; 227, 136, 137]/255;
shades{3} = [77, 175, 74; 88, 175, 85; 98, 175, 95; 107, 175, 105]/255;
shades{4} = [152, 78, 163; 155, 96, 163; 157, 114, 163; 159, 130, 163]/255;
shades{5} = [255, 127, 0; 255, 150, 55; 255, 171, 92; 255, 191, 128]/255;
shades{6} = [166, 86, 40; 167, 99, 60; 167, 112, 80; 166, 124, 100]/255;

shapes = {'o','diamond','square','>','<','*','^','v'};
nShapes = size(shapes,2);

% Affectation des couleurs aux groupes
if nGroup>nShapes*nColor
	for ng = 1:nShapes*nColor
		compGroup{ng}.colorGroup = shades{mod(ng-1,nColor)+1}(1,:);
	end
	palColorBlack = linspace(0,128,nGroup-nShapes*nColor)'*[1 1 1]/255; 
	for ng = nShapes*nColor+1:nGroup
		compGroup{ng}.colorGroup = palColorBlack(nGroup-ng+1,:);
	end
else
	for ng = 1:nGroup
		compGroup{ng}.colorGroup = shades{mod(ng-1,nColor)+1}(1,:);
	end
end

% Affectation des shapes aux groupes
for ng = 1:nGroup
	if ng<=nShapes*nColor
		compGroup{ng}.shapeGroup = shapes{ceil(ng/nColor)};
	else
		compGroup{ng}.shapeGroup = shapes{1};
	end
end

% Duplication des matrices facteurs pour passer en RGB
if boolRegroup==0
	A = y;
	for m = 1:M
		for r = 1:R
			A{m}(:,r) = A{m}(1:size(A{m},1),r)/max(A{m}(:,r));
		end
		A{m} = repmat(1-A{m}, 1, 1, 3);
	end
elseif boolRegroup==1
	A = cell(1,M);
	for m = 1:M
		A{m} = zeros(I,nGroup);
		for ng = 1:nGroup
			for r = 1:compGroup{ng}.sizeGroup
				A{m}(:,ng) = A{m}(:,ng)+lambda(compGroup{ng}.indGroup(r))*y{m}(:,compGroup{ng}.indGroup(r));
			end
			A{m}(:,ng) = A{m}(1:size(A{m},1),ng)/max(A{m}(:,ng));
		end
		A{m} = repmat(A{m}, 1, 1, 3);
	end
end

% Affichage des lambdas
figure(figNumber)
subplot(hRow*nRow+hLambda+hLambda+nRow+1,nCol,(nRow*hRow+hLambda+nRow+1)*nCol+[1:2*nCol])
hold on
for ng = 1:nGroup
	for r = 1:size(compGroup{ng}.indGroup,2)
		stem(compGroup{ng}.indGroup(r),100*lambda(compGroup{ng}.indGroup(r)),'filled','LineWidth',3,'Color',compGroup{ng}.colorGroup,'Marker',compGroup{ng}.shapeGroup,'MarkerSize',12)
	end
end
hold off
ax = gca;
ax.FontSize = 20;
if strcmpi(strLanguage,'EN')
	strLabelX = sprintf("Contribution of the %d components (in %%)",R);
elseif strcmpi(strLanguage,'FR')
	strLabelX = sprintf("Contribution des %d composantes (en %%)",R);
end
xlabel(strLabelX,'FontSize',30)
xlim([0 R+1])

indTicks = zeros(nGroup,1);
cellTicks = cell(1,nGroup);
for ng = 1:nGroup
    if compGroup{ng}.sizeGroup==1
        indTicks(ng) = inf;
    else
	    indTicks(ng) = compGroup{ng}.indGroup(ceil(0.5*compGroup{ng}.sizeGroup));
	    cellTicks{ng} = sprintf('%.1f%%',100*sum(lambda(compGroup{ng}.indGroup)));
    end
end
cellTicks = cellTicks(indTicks~=inf);
indTicks(indTicks==inf) = [];
xticks(indTicks)
xticklabels(cellTicks)

ax1 = cell(1,M);
% Affichage des Matrices facteurs
for row = 1:nRow
	for col = 1:nCol
		m = (row-1)*nCol+col;
		if m<M+1
			for ng = 1:nGroup
				colorGroupTemp = compGroup{ng}.colorGroup;
				if boolRegroup==0
					for r = 1:size(compGroup{ng}.indGroup,2)
						temp = y{m}(:,compGroup{ng}.indGroup(r));
						for i = 1:size(temp,2)
							temp(:,i) = temp(:,i)/max(temp(:,i));
						end
						A{m}(:,compGroup{ng}.indGroup(r),1) = 1-temp*(1-colorGroupTemp(1));
						A{m}(:,compGroup{ng}.indGroup(r),2) = 1-temp*(1-colorGroupTemp(2));
						A{m}(:,compGroup{ng}.indGroup(r),3) = 1-temp*(1-colorGroupTemp(3));
					end
				elseif boolRegroup==1
					temp = max(A{m}(:,ng,1));
					A{m}(:,ng,1) = 1-(A{m}(:,ng,1)/temp)*(1-colorGroupTemp(1));
					A{m}(:,ng,2) = 1-(A{m}(:,ng,2)/temp)*(1-colorGroupTemp(2));
					A{m}(:,ng,3) = 1-(A{m}(:,ng,3)/temp)*(1-colorGroupTemp(3));
				end
			end
			ax1{m} = subplot(hRow*nRow+hLambda+hLambda+nRow+1,nCol,(row-1)*nCol*(hRow+row-1)+(col-1)+[0:hRow-1]*nCol+1);
			imagesc(t{indVar(m)},1:R,permute(A{m},[2 1 3]));
			colormap(1-gray)
			ax = gca; ax.FontSize = 20;
			
			if xlimMan(1)<1e7
                if t{indVar(m)}(end)>1000
                    ecart = t{indVar(m)}(2)-t{indVar(m)}(1);
                    xlim([t{indVar(m)}(1)-ecart t{indVar(m)}(end)+ecart])
                else
				    xlim(xlimMan)
                end
			end
			
			ax = gca; ax.FontSize = 12;
			ytickangle(0);
			if strcmpi(strLanguage,'EN')
				if m==1, ylabel("Factor Matrices",'FontSize',30), end
			elseif strcmpi(strLanguage,'FR')
				if m==1, ylabel("Matrice Facteurs",'FontSize',40), end
			end
			title(strFluorescence{m},'FontSize',25)
		end
	end
end

linkaxes([ax1{1} ax1{2:end}],'y')

% Changement des couleurs du dendrogramme
subplot(hRow*nRow+hLambda+hLambda+nRow+1,nCol,(nRow*hRow+nRow)*nCol+[1:2*nCol])
rankToGroup = [];
for ng = 1:nGroup
	rankToGroup = [rankToGroup; ng*ones(length(compGroup{ng}.indGroup),1)];
end
zTorank = zeros(size(Z,1),1);

for z = 1:size(zTorank,1)-countBlack
	
	flag = 0;
	indGroup = z;
	while flag==0
		indGroup = min(Z(indGroup,1:2));
		if indGroup<R+1
			flag = 1;
		else
			indGroup = indGroup-R;
		end
	end
	zTorank(z) = rankToGroup(indGroup);
	if Z(z,3)<=colorThresh
		set(H(z),'Color',compGroup{zTorank(z)}.colorGroup);
	end
end

for z = 1:size(Z,1)
	if Z(z,3)>colorThresh
		set(H(z),'Color',[0 0 0]);
	end
end

for ng = 1:countBlack-1
	indNg = find(indPermute==nGroup-ng+1);
	plot(compGroup{indNg}.indGroup*[1 1], colorThresh*[0 1],'LineWidth',5,'Color',compGroup{indNg}.colorGroup);
end

xticks([1:R])
a = xticklabels();
xticks([1 5:5:(R-1) R])
xticklabels(a([1 5:5:(R-1) R]))

end

