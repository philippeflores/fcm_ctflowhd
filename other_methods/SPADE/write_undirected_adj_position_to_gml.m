function [node_names] = write_undirected_adj_position_to_gml(filename, adj, node_positions, node_names)
% write_undirected_adj_position_to_gml(filename, adj, positions, node_names)
% 
% the network must be undirected

if ~exist('node_names')
    node_names= cell(1,size(adj,1));
    for i=1:length(node_names)
        node_names{i} = ['node-',num2str(i)];
    end
end


fid = fopen(filename,'w');
fprintf(fid,'Creator\t"Cytoscape"\n');
fprintf(fid,'Version\t1.0\n');
fprintf(fid,'graph [\n');

all_dist = pdist(node_positions');
node_positions = round(node_positions./min(all_dist(all_dist~=0))*20);
node_positions(2,:) = -node_positions(2,:);

node_to_write=[];
for i=1:size(adj,1)
    str_of_this_node = ['\tnode\t[\n', ...
                        ['\t\troot_index\t',num2str(-size(adj,1)+i-1),'\n'], ...
                        ['\t\tid\t', num2str(-size(adj,1)+i-1),'\n'], ...
                        '\t\tgraphics\t[\n', ...
                        '\t\t\tx\t', num2str(node_positions(1,i)),'\n', ...
                        '\t\t\ty\t', num2str(node_positions(2,i)),'\n', ...
                        '\t\t\tw\t40\n', ...
                        '\t\t\th\t40\n', ...
                        '\t\t\tfill\t"#ff9999"\n',...
                        '\t\t\ttype\t"ellipse"\n',...
                        '\t\t\toutline\t"#666666"\n',...
                        '\t\t\toutline_width\t1.5\n',...
                        '\t\t]\n',...
                        '\t\tlabel\t"',node_names{i},'"\n',...
                        '\t]\n',...
                        ] ;
    fprintf(fid,str_of_this_node);
end




ind = find(triu(adj,1)'==1);
edge_end = mod(ind, size(adj,1)); edge_end(edge_end==0) = size(adj,1);
edge_start = (ind - edge_end)/size(adj,1)+1;
for i=1:length(edge_start)
    str_of_this_edge = ['\tedge\t[\n', ...
                        ['\t\troot_index\t',num2str(-length(edge_start)+i-1),'\n'], ...
                        ['\t\ttarget\t', num2str(-size(adj,1)+edge_end(i)-1),'\n'], ...
                        ['\t\tsource\t', num2str(-size(adj,1)+edge_start(i)-1),'\n'], ...
                        '\t\tgraphics\t[\n', ...
                        '\t\t\twidth\t1.5\n', ...
                        '\t\t\tfill\t"#0000ff"\n',...
                        '\t\t\ttype\t"line"\n',...
                        '\t\t\tLine\t[\n',...
                        '\t\t\t]\n',...
                        '\t\t\tsource_arrow\t0\n',...
                        '\t\t\ttarget_arrow\t0\n',...
                        '\t\t]\n',...
                        '\t\tlabel\t"edge"\n',...
                        '\t]\n',...
                        ] ;
    fprintf(fid,str_of_this_edge);
end
fprintf(fid,']\n');
fclose(fid);



% left_columns = cell(sum(sum(triu(adj))),3);
% counter = 1;
% for i=1:size(adj,1)
%     for j=i:size(adj,1)
%         if adj(i,j)==1
%             left_columns(counter,1) = node_names(i);
%             left_columns(counter,2) = {'edge'};
%             left_columns(counter,3) = node_names(j);
%             counter = counter+1;
%         end
%     end
% end
% write_to_txt(filename, [], left_columns, [], char(9));
