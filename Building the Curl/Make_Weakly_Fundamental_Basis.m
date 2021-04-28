function [C,chords,chord_order,search_cost] = make_simple_basis(adjacency,tree_mode)

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

%% make endpoints to edge index mapping
edge_index = sparse(endpoints(:,1),endpoints(:,2),(1:E),V,V) - sparse(endpoints(:,2),endpoints(:,1),(1:E),V,V);

%% generate vector of all zeros since we'll use it repeatedly
z = zeros(V,1);

if connected == 1
    %% initialize search network, stored in search adjacency
    search_adjacency = cell(V,1);
    for k = 1:V-1
        edge = tree(k);
        nodes = endpoints(edge,:);
        search_adjacency{nodes(1)} = [search_adjacency{nodes(1)},nodes(2)];
        search_adjacency{nodes(2)} = [search_adjacency{nodes(2)},nodes(1)];
    end
    
    %% build ancestor lists to get distance of each endpoint to root
    [~,order] = sort(node_indices,'ascend');
    ancestor_edges = cell(V,1);
    distance = z;
    for j = 2:V
        node = order(j);
        parent = parents(node);
        parent_edge = parent_edges(node);
        ancestor_edges{node} = [ancestor_edges{parent},parent_edge];
        distance(node) = length(ancestor_edges{node});
    end
    
    %% sort chords by distance of further endpoint (chords should already be ordered this way if use a breadth first search)
    chord_distances = zeros(L,1);
    for k = 1:L
        chord = chords(k);
        nodes = endpoints(chord,:);
        chord_distances(k) = mean(distance(nodes)); % all chords either are in the same generation, or pass between adjacent generations
    end
    [~,chord_order] = sort(chord_distances,'ascend');
    
    %% preallocate list of edges in loops
    edges_list = []; % will store all edges in loops
    directions_list = []; % will store direction the edges are crossed
    
    %% loop over chords
    search_cost = 0;
    for k = 1:L
        %% find chord an endpoints
        chord = chords(chord_order(k));
        
        nodes = endpoints(chord,:);
        root_1 = nodes(1);
        root_2 = nodes(2);
        
        %% initialize search
        leaves_1 = root_1;
        leaves_2 = root_2;
        
        found_1 = z; % tracks which nodes have been reached by the search from root 1
        found_2 = z; % tracks which nodes have been reached by the search from root 2
        found_1(leaves_1) = 1;
        found_2(leaves_2) = 1;
        
        parents_1 = z;
        parents_2 = z;
        
        intersection_nodes = []; % will store all nodes where the search trees intersect
        
        stop = 0;
        while stop == 0
            %% get parents
            parent_1 = leaves_1(1);
            parent_2 = leaves_2(1);
            
            leaves_1 = leaves_1(2:end);
            leaves_2 = leaves_2(2:end);
            
            %% get neighbors
            neighbors_1 = search_adjacency{parent_1};
            neighbors_2 = search_adjacency{parent_2};
            
            %% loop over neighbors
            for j = 1:length(neighbors_1)
                neighbor = neighbors_1(j);
                
                %% update search cost
                search_cost = search_cost + 1;
                
                %% check if we've seen it before
                if found_1(neighbor) == 0
                    %% it is new
                    found_1(neighbor) = 1;
                    leaves_1 = [leaves_1,neighbor];
                    
                    parents_1(neighbor) = parent_1;
                    
                    %% check if found by other search
                    if found_2(neighbor) == 1
                        stop = 1;
                        intersection_nodes = [intersection_nodes,neighbor];
                    end
                end
            end
            
            for j = 1:length(neighbors_2)
                neighbor = neighbors_2(j);
                
                %% check if we've seen it before
                if found_2(neighbor) == 0
                    %% it is new
                    found_2(neighbor) = 2;
                    leaves_2 = [leaves_2,neighbor];
                    
                    parents_2(neighbor) = parent_2;
                    
                    %% check if found by other search
                    if found_1(neighbor) == 1
                        stop = 1;
                        intersection_nodes = [intersection_nodes,neighbor];
                    end
                end
            end
            
        end
        
        %% pick a random intersection to start from
        intersection = intersection_nodes(randperm(length(intersection_nodes),1));
        
        %% trace back through each tree
        stop = 0;
        edges_1 = [];
        directions_1 = [];
        
        start = intersection;
        if start ~= root_1
            while stop == 0
                finish = parents_1(start);
                edge = edge_index(start,finish);
                directions_1 = [directions_1,sign(edge)];
                edges_1 = [edges_1,abs(edge)];
                
                if finish == root_1
                    stop = 1;
                else
                    start = finish;
                end
            end
        end
        
        
        
        stop = 0;
        edges_2 = [];
        directions_2 = [];
        
        start = intersection;
        if start ~= root_2
            while stop == 0
                finish = parents_2(start);
                
                edge = edge_index(start,finish);
                directions_2 = [directions_2,sign(edge)];
                edges_2 = [edges_2,abs(edge)];
                
                if finish == root_2
                    stop = 1;
                else
                    start = finish;
                end
            end
        end
        
        %% build cycle
        edges_in_cycle = [chord,edges_2,edges_1];
        directions = [1,directions_2,-directions_1];
        
        %% add chord to search network
        search_adjacency{root_1} = [search_adjacency{root_1},root_2];
        search_adjacency{root_2} = [search_adjacency{root_2},root_1];
        
        %% add to list
        edges_list = [edges_list, edges_in_cycle + (k - 1)*E];
        directions_list = [directions_list, directions];
        
    end
    
    %% build curl
    C = reshape(sparse(edges_list,1,directions_list,L*E,1),[E,L])';

else
    C = nan;
end