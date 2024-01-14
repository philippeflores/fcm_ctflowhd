function plotFigSpade(distR,y,lambda,strFluorescence,compGroup,varargin)

strScreen = 's';
% Initialisation des variables
nGroup = size(compGroup,2);
R = size(distR,1);
M = size(y,2);
nMax = 6;
I = size(y{1},1);

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRSCREEN'
                strScreen = varargin{i+1};
            case 'NMAX'
                nMax = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end



k = floor((M>=(nMax-1))*1+1);
if M>10
	L = ceil((M-(nMax*2-4))/nMax)+2;
else
	L = ceil((M+k*k)/nMax);
end
C = min(nMax,k+ceil(M/L));
indK = zeros(k*k,1); for ind = 1:k*k, modind = mod(ind,k); modind(modind==0)=k; indK(ind) = (floor((ind-1)/k)*C)+modind; end
indReste = 1:(L*C); indReste(indK)=[];

[adj,~,~] = mst(distR);
node_positions = arch_layout(adj);

if strcmpi(strScreen,'s')
	vScreen = [0 0 1 1];
elseif strcmpi(strScreen,'b')
	vScreen = [-1.43, 0, 1.43,1.33];
else
	vScreen = [0 0 1 1];
end

figure
set(gcf,'Units','Normalized','OuterPosition',vScreen)

% Affichage des résultats
ax = subplot(L,C,indK);
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
		a = plot(node_positions(1,structGroup.indGroup(orderPlot(r))),node_positions(2,structGroup.indGroup(orderPlot(r))),'MarkerFaceColor',structGroup.colorGroup,'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*25*R*k*k/(L*C),'MarkerEdgeColor','k','Marker',compGroup{ng}.shapeGroup);
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
	xticks([]), yticks([])
end

axExp = cell(1,M);
cmapJet = colormap([50,136,189; 102,194,165; 171,221,164; 230,245,152; 255,255,191; 254,224,139; 253,174,97; 244,109,67; 213,62,79]/255);
for m = 1:M
	axExp{m} = subplot(L,C,indReste(m));
	hold on
	for r = 1:R
		indConnex = find(adj(r,:));
		for ind = 1:size(indConnex,2)
			plot([node_positions(1,indConnex(ind)) node_positions(1,r)],[node_positions(2,indConnex(ind)) node_positions(2,r)],'-k','LineWidth',2)
		end
	end
	posEsp = round(sum(y{m}.*(repmat(1:I,R,1)')));
% 	posEsp = (posEsp-min(posEsp))+1;
% 	posEsp = round(posEsp*I/max(posEsp));
	for ng = 1:nGroup
		structGroup = compGroup{ng};
		[~,orderPlot] = sort(lambda(structGroup.indGroup),'descend');
		for r = 1:size(structGroup.indGroup,2)
			a = plot(node_positions(1,structGroup.indGroup(orderPlot(r))),node_positions(2,structGroup.indGroup(orderPlot(r))),'MarkerFaceColor',cmapJet(ceil(posEsp(structGroup.indGroup(orderPlot(r)))*size(cmapJet,1)/I),:),'MarkerSize',(lambda(structGroup.indGroup(orderPlot(r)))+1/R)*25*k*R/(L*C),'MarkerEdgeColor','k','Marker',compGroup{ng}.shapeGroup);
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
% 	cb = colorbar('Location','southoutside','Ticks',[]);
	axExp{m}.FontSize = 25;
    clim([-0.5 4.5])
    colorbar
	xticks([]), yticks([])
	title(sprintf("%s",strFluorescence{m}),'FontSize',35)
	
end
linkaxes([ax axExp{:}],'xy')

end

function [adj,adj2, cost_value] = mst(X,working_mode,exclude_adj)

% Minimal or Minimum Spanning Tree based on Euclidian distances
% MST in short: use (X)(n x p) to form (n-1) lines to connect (n) objects in the shortest possible way in the (p)
% dimensional variable-space, under the condition 'no closed loops allowed'.
% working_mode : 'euclidean' (default)
%                'corr'
%                'abs_corr'
% 
% out:Xmst (objects-1 x 2) link set between 'objects' indexed as rows in X
%     adj adjacency matrix


cost_value = 0;

if ~exist('working_mode','var')
    working_mode = 'euclidean';
end
if isempty(intersect({'euclidean','corr','abs_corr'}, working_mode))
    working_mode = 'euclidean';
end
if ~isempty(intersect({'corr','abs_corr'}, working_mode))
    X = per_gene_normalization(X); % this makes computing correlation easier
    X = X./norm(X(1,:));
end

if ~exist('exclude_adj','var') || isempty(exclude_adj) 
    exclude_adj = sparse(size(X,1),size(X,1));
end

nX = size(X,1);

components = []; active_components =[];
adj = sparse(size(X,1),size(X,1)); 
adj2 = sparse(size(X,1),size(X,1)); 
count = 0; %fprintf('constructing a total of %d MST edges ... %6d', size(X,1)-1,count);
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
        count = count + 1; %fprintf('\b\b\b\b\b\b%6d', count);
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
%     sum(active_components)
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
    count = count + 1; %fprintf('\b\b\b\b\b\b%6d', count);
end
% fprintf('\n');
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
% dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-corr; 
return
end

function dist = comp_dist_abs_corr(X,ind1,ind2)
% dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-abs(corr); 
return
end

function node_positions = arch_layout(tree_adj)
warning off;
% tree_adj is undirected
adj = triu(tree_adj,1); adj = adj + adj';

% % get shortest hop
% p = adj+eye(size(adj));
% sh = p; hop=2; new_links = sh;
% while(sum(sum(sh==0))~=0), 
%     hop
%     new_links = (p*new_links~=0 & sh==0);
%     sh = sh + new_links*hop; hop = hop+1; 
% end % this line replaces the shortest hop calculation
% shortest_hop=sh;

shortest_hop = tree_shortest_hop(adj);
shortest_hop = shortest_hop - diag(diag(shortest_hop));


% find the backbone, which is the longest path
[ind_j,ind_i]  = find_matrix_top_element(shortest_hop);
ind_i = ind_i(1); ind_j = ind_j(1);
back_bones = find(shortest_hop(ind_i,:)+shortest_hop(ind_j,:)==shortest_hop(ind_i,ind_j));
[~,I] = sort(shortest_hop(ind_i,back_bones));
back_bones = back_bones(I);

% find all the side_chains
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
        % find the side chain starting from root_node to first_neighbor(i) to as far as it can go
        % 1 find members on this side chain and the subside chains
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
% determin the nodes location of the backbones
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
% DrawTree(adj(position_assigned_flag==1,position_assigned_flag==1), coeff(:,position_assigned_flag==1)); drawnow
% fprintf('Drawing a total of %d nodes ... %6d', size(tree_adj,1),sum(position_assigned_flag==1));

% determin the order of subbackbones to draw
centernees = [];
for k=1:length(side_chains)
    if ismember(side_chains{k}(1),back_bones)==0, break; end
    centernees(k) = abs(shortest_hop(side_chains{k}(1),back_bones(1)) - shortest_hop(side_chains{k}(1),back_bones(end)));
end
[~,I] = sort(centernees);
% determin the nodes location of each node in the side chains
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
%                 repel_force = sum(repel_vector./repmat(sqrt(sum(repel_vector.^2,1)).^5,2,1).*repmat(sum(adj(:,position_assigned_flag==1)~=0),2,1),2);
                
                cos_alfa = (repel_force'*[cos(theta(n));sin(theta(n))])/norm(repel_force);
                sin_alfa = sqrt(1-cos_alfa^2);
                potential_position_force(m,n,1) = norm(repel_force)*sin_alfa; % force perpendicular to string
                potential_position_force(m,n,2) = norm(repel_force)*cos_alfa; % force on string
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
%         DrawTree(adj(position_assigned_flag==1,position_assigned_flag==1), coeff(:,position_assigned_flag==1)); drawnow
%         fprintf('\b\b\b\b\b\b%6d', sum(position_assigned_flag==1));
    end
end
node_positions = node_positions - repmat(median(node_positions,2),1,size(node_positions,2));
node_positions(1,:) = node_positions(1,:)/max(abs(node_positions(1,:)))*50;
node_positions(2,:) = node_positions(2,:)/max(abs(node_positions(2,:)))*50;
% fprintf('\n');

return
end



% function DrawTree(adj, node_positions, node_sizes)
% 
% coeff = node_positions;
% 
% hold off; plot(0); hold on;
% pairs = find_matrix_big_element(triu(adj,1),1);
% for k=1:size(pairs,1), line(coeff(1,pairs(k,:)),coeff(2,pairs(k,:)),'color','g'); end
% % % draw labels
% % for k=1:size(coeff,2), text(coeff(1,k)+2,coeff(2,k),num2str(k),'FontSize',7); end
% 
% if exist('node_sizes','var')
%     log_node_size = log10(node_size); log_node_size(log_node_size<=0) = min(log_node_size>0);
%     draw_node_size = round(log_node_size/median(log_node_size)*8);
%     draw_node_size(draw_node_size==0) = 1;
% else
%     draw_node_size = ones(1,size(coeff,2))*8;
% end
% 
% for k=1:size(coeff,2) 
%     plot(coeff(1,k),coeff(2,k),'o','markersize',draw_node_size(k),'markerfacecolor',[0 0 0],'color',[0 0 0]); 
% end
% hold off;
% axis(reshape([-max(abs(coeff)');+max(abs(coeff)')],1,4)*1.1);
% 
% end


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
% shortest_hop = tree_shortest_hop(adj_matrix)

adj_matrix = (abs(adj_matrix) + abs(adj_matrix') + eye(size(adj_matrix))) >0;
nNodes=size(adj_matrix,1);
shortest_hop = zeros(nNodes);
e = zeros(nNodes,1);
e(1)=1;

while sum(e==0)~=0  % e vector serves as a flag, see whether all nodes are "included/esamined"
    ind = find(sum(adj_matrix(e==1,:),1)'~=0 & e==0); % ind of all the nodes that are connected to one of the elements in e
    ind_new = ind(1); % take one of them
    ind_exist = find(adj_matrix(:,ind(1))==1 & e==1); % the one node in e that the new node connects to
    shortest_hop(ind_new,e==1) = shortest_hop(ind_exist,e==1)+1;
    shortest_hop(e==1,ind_new) = shortest_hop(e==1,ind_exist)+1;
    e(ind_new)=1;
end

shortest_hop = shortest_hop + eye(size(shortest_hop));
end







