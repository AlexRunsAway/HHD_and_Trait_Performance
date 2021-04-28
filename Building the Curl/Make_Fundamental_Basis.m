function C = make_fundamental_basis(adjacency,tree_mode)

%% build spanning tree and find chords
if strcmp(tree_mode,'depth')
    output = depth_first_search(adjacency);
elseif strcmp(tree_mode,'degree')
    output = degree_first_search(adjacency);
else
    output = breadth_first_search(adjacency);
end

node_indices = output.node_index; % order nodes were found in search
endpoints = output.endpoints; % endpoints of each edge
tree = output.tree; % edges in tree
parent_edges = output.parent_edges;
parents = output.parents;
chords = output.chords; % chords
connected = output.connected; % logical flag for whether the graph is connected

%% get dimensions
V = length(node_indices);
[E,~] = size(endpoints);
L = E - (V - 1);

if connected == 1
    %% build ancestor lists
    [~,order] = sort(node_indices,'ascend');
    ancestor_edges = cell(V,1);
    for j = 2:V
        node = order(j);
        parent = parents(node);
        parent_edge = parent_edges(node);
        ancestor_edges{node} = [ancestor_edges{parent},parent_edge];
    end
    
    %% construct curl
    paths.root_to_start = cell(1,L);
    paths.finish_to_root = cell(1,L);
    for l = 1:L
        
        chord = chords(l);
        nodes = endpoints(chord,:);
        start = nodes(1);
        finish = nodes(2);
        
        paths.root_to_start{l} = ancestor_edges{start} + (l - 1)*E;
        paths.finish_to_root{l} = ancestor_edges{finish} + (l - 1)*E;
    end
    
    Paths.root_to_start = cell2mat(paths.root_to_start);
    Paths.finish_to_root = cell2mat(paths.finish_to_root);
    
    Out_vector = sparse(Paths.root_to_start,1,1,L*E,1);
    Back_vector = sparse(Paths.finish_to_root,1,1,L*E,1);
    
    Out = reshape(Out_vector,[E,L])';
    Back = reshape(Back_vector,[E,L])';
    
    % Construct curl
    C = Out - Back + sparse((1:L),chords,1,L,E);

else
    C = nan;
end