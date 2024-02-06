function plot_MST(y,lambda,strLabel,compGroup,t,indVar,varargin)

strScreen = 'halfR';
colMax = 2;
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
            case 'STRCENTROID'
                strCentroid = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

colMax = colMax+2;
M = size(y,2);
I = size(y{1},1);
nGroup = size(compGroup,2);
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

k = 2;
if M>10
	nRow = ceil((M-(colMax*2-4))/colMax)+2;
else
	nRow = ceil((M+k*k)/colMax);
end
nCol = min(colMax,k+ceil(M/nRow));
indK = zeros(k*k,1); for ind = 1:k*k, modind = mod(ind,k); modind(modind==0)=k; indK(ind) = (floor((ind-1)/k)*nCol)+modind; end
indRem = 1:(nRow*nCol); indRem(indK)=[];

[adj,~,~] = mst(matFeature);
node_positions = arch_layout(adj);

figure
bigScreen(strScreen)

ax = subplot(nRow,nCol,indK);
hold on 
for r = 1:R
	indConnex = find(adj(r,:));
	for ind = 1:size(indConnex,2)
		plot([node_positions(1,indConnex(ind)) node_positions(1,r)],[node_positions(2,indConnex(ind)) node_positions(2,r)],'-k','LineWidth',2)
	end
end

for ng = 1:nGroup
	structGroup = compGroup{ng};
	[~,orderPlot] = sort(lambda(structGroup.indGroup),'descend');
	for r = 1:size(structGroup.indGroup,2)
		a = plot(node_positions(1,structGroup.indGroup(orderPlot(r))),node_positions(2,structGroup.indGroup(orderPlot(r))),'MarkerFaceColor',structGroup.colorGroup,'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*25*R*k*k/(nRow*nCol),'MarkerEdgeColor','k','Marker',compGroup{ng}.shapeGroup);
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
	title("Hierarchical clustering",'FontSize',35)
    ax.FontSize = 15;
end

axExp = cell(1,M);
cmapJet = colormap([50,136,189; 102,194,165; 171,221,164; 230,245,152; 255,255,191; 254,224,139; 253,174,97; 244,109,67; 213,62,79]/255);
for m = 1:M
	axExp{m} = subplot(nRow,nCol,indRem(m));
	hold on
	for r = 1:R
		indConnex = find(adj(r,:));
		for ind = 1:size(indConnex,2)
			plot([node_positions(1,indConnex(ind)) node_positions(1,r)],[node_positions(2,indConnex(ind)) node_positions(2,r)],'-k','LineWidth',2)
		end
	end
	posEsp = round(sum(y{m}.*(repmat(1:I,R,1)')));
	for ng = 1:nGroup
		structGroup = compGroup{ng};
		[~,orderPlot] = sort(lambda(structGroup.indGroup),'descend');
		for r = 1:size(structGroup.indGroup,2)
			a = plot(node_positions(1,structGroup.indGroup(orderPlot(r))),node_positions(2,structGroup.indGroup(orderPlot(r))),'MarkerFaceColor',cmapJet(ceil(posEsp(structGroup.indGroup(orderPlot(r)))*size(cmapJet,1)/I),:),'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*25*k*R/(nRow*nCol),'MarkerEdgeColor','k','Marker',compGroup{ng}.shapeGroup);
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
	axExp{m}.FontSize = 15;
    clim([-0.5 4.5])
    cb = colorbar;
    cb.Label.Position(1) = 4;
    if t{indVar(m)}(end)<100
        ylabel(cb,'Fluorescence','FontSize',18,'Rotation',270)
    end
	title(sprintf("%s",strLabel{m}),'FontSize',25)
	
end
linkaxes([ax axExp{:}],'xy')

end

function [adj,adj2, cost_value] = mst(X,working_mode,exclude_adj)

cost_value = 0;

if ~exist('working_mode','var')
    working_mode = 'euclidean';
end
if isempty(intersect({'euclidean','corr','abs_corr'}, working_mode))
    working_mode = 'euclidean';
end
if ~isempty(intersect({'corr','abs_corr'}, working_mode))
    X = per_gene_normalization(X);
    X = X./norm(X(1,:));
end

if ~exist('exclude_adj','var') || isempty(exclude_adj) 
    exclude_adj = sparse(size(X,1),size(X,1));
end

nX = size(X,1);

components = []; active_components =[];
adj = sparse(size(X,1),size(X,1)); 
adj2 = sparse(size(X,1),size(X,1)); 
count = 0;
for i=1:nX
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(X,i,1:nX); 
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(X,i,1:nX);
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(X,i,1:nX); 
    end
    dist(i) = max(dist)+1; 
    dist = dist + exclude_adj(i,:).*(max(dist)+1);    
    dist = full(dist);
    [Dmin,Dwin] = min(dist);
    Xmst(i,:) = [i Dwin];
    if adj(i,Dwin)==0 && adj(Dwin,i)==0
        adj(i,Dwin)=1;adj(Dwin,i)=1;
        adj2(i,Dwin)=Dmin;adj2(Dwin,i)=Dmin;
        cost_value = cost_value + Dmin;
        count = count + 1;
    end
    if isempty(components)
        components = sparse(zeros(size(X,1),1)); components([i,Dwin],1) = 1; active_components=1;
    else
        [existing_comp1] = find(components(i,:)==1 & active_components==1);
        [existing_comp2] = find(components(Dwin,:)==1 & active_components==1);
        if isempty(existing_comp1) && isempty(existing_comp2)
            components = [components,zeros(size(X,1),1)]; components([i,Dwin],end) = 1; active_components = [active_components,1];
        elseif ~isempty(existing_comp1) && isempty(existing_comp2)
            components([i,Dwin],existing_comp1)=1;
        elseif isempty(existing_comp1) && ~isempty(existing_comp2)
            components([i,Dwin],existing_comp2)=1;
        elseif ~isempty(existing_comp1) && ~isempty(existing_comp2) && existing_comp1~=existing_comp2
            components = [components, components(:,existing_comp1)+components(:,existing_comp2)];
            active_components = [active_components,1];
            active_components([existing_comp1, existing_comp2])=0;
        end
    end
end
 
while sum(active_components)>1
    components_sizes = sum(components); components_sizes(active_components==0) = max(components_sizes+1);
    [~, existing_comp1] = min(components_sizes);
    ind1 = find(components(:,existing_comp1)==1); ind1 = ind1(:)';
    ind2 = setdiff(1:size(components,1),ind1); ind2 = ind2(:)';
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(X,ind1,ind2); 
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(X,ind1,ind2);
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(X,ind1,ind2);
    end
    dist = dist + exclude_adj(ind1,ind2).*(max(max(dist))+1); dist = full(dist);
    [Dmin,ind] = min(reshape(dist,length(ind1)*length(ind2),1));
    j = ceil(ind/length(ind1));
    i = ind - (j-1)*length(ind1);
    Xmst = [Xmst; [ind1(i),ind2(j)]];
    adj(ind1(i),ind2(j))=1; adj(ind2(j),ind1(i))=1; 
    adj2(ind1(i),ind2(j))=Dmin; adj2(ind2(j),ind1(i))=Dmin; 
    cost_value = cost_value + Dmin;
    [existing_comp2] = find(components(ind2(j),:)==1 & active_components==1);
    components(:,existing_comp1) = components(:,existing_comp1) + components(:,existing_comp2);
    active_components(existing_comp2)=0;
    count = count + 1;
end
return

end


function dist = comp_dist_euclidean(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
for i=1:length(ind1)
    dist(i,:) = sqrt(sum((repmat(X(ind1(i),:),length(ind2),1) - X(ind2,:)).^2,2)); 
end
return
end

function dist = comp_dist_corr(X,ind1,ind2)
corr = X(ind1,:)*X(ind2,:)';
dist = 1-corr; 
return
end

function dist = comp_dist_abs_corr(X,ind1,ind2)
corr = X(ind1,:)*X(ind2,:)';
dist = 1-abs(corr); 
return
end

function node_positions = arch_layout(tree_adj)
warning off;
adj = triu(tree_adj,1); adj = adj + adj';

shortest_hop = tree_shortest_hop(adj);
shortest_hop = shortest_hop - diag(diag(shortest_hop));

[ind_j,ind_i]  = find_matrix_top_element(shortest_hop);
ind_i = ind_i(1); ind_j = ind_j(1);
back_bones = find(shortest_hop(ind_i,:)+shortest_hop(ind_j,:)==shortest_hop(ind_i,ind_j));
[~,I] = sort(shortest_hop(ind_i,back_bones));
back_bones = back_bones(I);

side_chains = cell(0);
side_chain_roots = back_bones;
counter=1;
while counter<=length(side_chain_roots)
    root_node = side_chain_roots(counter);
    first_neighbors = setdiff(find(adj(root_node,:)~=0),side_chain_roots);
    if isempty(first_neighbors)
        counter = counter + 1;
        continue;
    end
    for i=1:length(first_neighbors)
        subtree_nodes_through_this_neighbor = find(shortest_hop(root_node,:)>shortest_hop(first_neighbors(i),:));
        [~,I] = max(shortest_hop(root_node,subtree_nodes_through_this_neighbor));
        end_node = subtree_nodes_through_this_neighbor(I);
        sub_backbone = find(shortest_hop(root_node,:)+shortest_hop(end_node,:)==shortest_hop(root_node,end_node));
        [~,I] = sort(shortest_hop(root_node,sub_backbone));
        sub_backbone = sub_backbone(I);
        
        side_chain_roots = [side_chain_roots, sub_backbone(2:end)];
        side_chains = [side_chains;{sub_backbone}];
    end
end


node_positions = zeros(2,size(adj,1));
position_assigned_flag = zeros(1,size(adj,1));
backbone_node_angles = 1:length(back_bones);
backbone_node_angles = backbone_node_angles - mean(backbone_node_angles);
backbone_node_angles = backbone_node_angles./max(backbone_node_angles).*(pi/90*25);
node_positions(1,back_bones) = sin(backbone_node_angles);
node_positions(2,back_bones) = -cos(backbone_node_angles);
node_positions = node_positions./norm(node_positions(:,back_bones(1))-node_positions(:,back_bones(2)));
if length(back_bones)>500
    node_positions = node_positions.*(500/length(back_bones));
end
position_assigned_flag(back_bones)=1;


coeff = node_positions - repmat(median(node_positions,2),1,size(node_positions,2));
coeff(1,:) = coeff(1,:)/max(abs(coeff(1,:)))*50;
coeff(2,:) = coeff(2,:)/max(abs(coeff(2,:)))*50;

centernees = [];
for k=1:length(side_chains)
    if ismember(side_chains{k}(1),back_bones)==0, break; end
    centernees(k) = abs(shortest_hop(side_chains{k}(1),back_bones(1)) - shortest_hop(side_chains{k}(1),back_bones(end)));
end
[~,I] = sort(centernees);
for i=[I,k:length(side_chains)]
    for j=1:length(side_chains{i})
        if position_assigned_flag(side_chains{i}(j))==1
            continue;
        end
        new_node = side_chains{i}(j);
        attaching_node = find(adj(new_node,:)~=0 & position_assigned_flag==1);
        r = 0.3:0.1:0.9;
        theta = 0:2*pi/90:2*pi;
        potential_position_force = zeros(length(r),length(theta),2);
        for m = 1:length(r)
            for n = 1:length(theta)
                potential_position = node_positions(:,attaching_node) + r(m)*[cos(theta(n));sin(theta(n))];
                repel_vector = repmat(potential_position,1,sum(position_assigned_flag==1 & shortest_hop(new_node,:)<200)) - node_positions(:,position_assigned_flag==1 & shortest_hop(new_node,:)<200);
                repel_force = sum(repel_vector./repmat(sqrt(sum(repel_vector.^2,1)).^5,2,1),2);
                
                cos_alfa = (repel_force'*[cos(theta(n));sin(theta(n))])/norm(repel_force);
                sin_alfa = sqrt(1-cos_alfa^2);
                potential_position_force(m,n,1) = norm(repel_force)*sin_alfa; 
                potential_position_force(m,n,2) = norm(repel_force)*cos_alfa; 
            end
        end
        best_for_each_layer = []; best_string_force_each_layer=[];
        for m = 1:length(r)
            n_s = find(potential_position_force(m,:,2)>=0);
            if isempty(n_s), continue; end
            [~,I] = min(potential_position_force(m,n_s,1));
            best_for_each_layer = [best_for_each_layer;[m,n_s(I)]];
            best_string_force_each_layer = [best_string_force_each_layer; potential_position_force(m,n_s(I),2)];
        end
        [~,I] = min(best_string_force_each_layer);
        best_r = r(best_for_each_layer(I,1));
        best_theta = theta(best_for_each_layer(I,2));
        
        node_positions(:,new_node) = node_positions(:,attaching_node) + best_r*[cos(best_theta);sin(best_theta)];
        position_assigned_flag(side_chains{i}(j))=1;
        
        coeff = node_positions - repmat(median(node_positions,2),1,size(node_positions,2));
        coeff(1,:) = coeff(1,:)/max(abs(coeff(1,:)))*50;
        coeff(2,:) = coeff(2,:)/max(abs(coeff(2,:)))*50;
    end
end
node_positions = node_positions - repmat(median(node_positions,2),1,size(node_positions,2));
node_positions(1,:) = node_positions(1,:)/max(abs(node_positions(1,:)))*50;
node_positions(2,:) = node_positions(2,:)/max(abs(node_positions(2,:)))*50;

return
end




function [ind_i,ind_j] = find_matrix_top_element(adj)

ind = find(adj==max(max(adj)));
dim = size(adj,1);
for i=1:length(ind)
    ind_j(i) = ceil(ind(i)/dim);
    ind_i(i) = ind(i) - (ind_j(i)-1)*dim;
end
return
end

function shortest_hop = tree_shortest_hop(adj_matrix)

adj_matrix = (abs(adj_matrix) + abs(adj_matrix') + eye(size(adj_matrix))) >0;
nNodes=size(adj_matrix,1);
shortest_hop = zeros(nNodes);
e = zeros(nNodes,1);
e(1)=1;

while sum(e==0)~=0
    ind = find(sum(adj_matrix(e==1,:),1)'~=0 & e==0); 
    ind_new = ind(1); 
    ind_exist = find(adj_matrix(:,ind(1))==1 & e==1);
    shortest_hop(ind_new,e==1) = shortest_hop(ind_exist,e==1)+1;
    shortest_hop(e==1,ind_new) = shortest_hop(e==1,ind_exist)+1;
    e(ind_new)=1;
end

shortest_hop = shortest_hop + eye(size(shortest_hop));
end







