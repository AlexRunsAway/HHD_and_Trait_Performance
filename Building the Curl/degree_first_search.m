function output = degree_first_search(adjacency)
%% inputs
% 1. adjacency is an adjacency structure

%% outputs
% 0. node indices in tree (order nodes discovered)
% 1. edge to endpoints mapping
% 2. list of edges in tree
% 3. list of parent edges
% 4. list of chords
% 5. connected, if connected = 0 then the graph is not connected

%% get number of nodes
V = length(adjacency);

%% find degree of each node
degree = nan(V,1);
for j = 1:V
    degree(j) = length(adjacency{j});
end

%% sort the nodes by their degree
[~,degree_order] = sort(degree,'descend');

%% find number of edges and loops
E = sum(degree)/2;
L = E - (V - 1);

%% preallocate
endpoints = nan(E,2); % will store the endpoints of each edge

node_index = nan(V,1); % will store whether or not a node has been found by the search procedure yet
parent_edges = nan(V,1); % will store the parent edge for each node
parents = nan(V,1); % will store the parent of each node

tree = nan(V-1,1); % will store the edges in the tree
chords = nan(L,1); % will store the chords

%% initialize search
edge_index = 0;
chords_found = 0;
searched_from = nan(V,1); % equals nan if we haven't searched from this node, 1 otherwise
tree_size = 1;
leaves = degree_order(1); % nodes we could still search from


%% search
stop = 0;
while stop == 0
    
     %% update current node (set to first leaf), remove from list of leaves
    parent = leaves(1);
    node_index(parent) = tree_size;
    searched_from(parent) = 1;
    
    leaves = leaves(2:end); 
    
    %% find neighbors
    neighbors = adjacency{parent};
    
    %% loop over neighbors, add to leaves if not found yet
    for k = 1:length(neighbors)
        neighbor = neighbors(k);
        
        if neighbor ~= parents(parent) % don't backtrack
            if isnan(node_index(neighbor)) % check if already found
                %% a new node has been found
                % store parent of child
                child = neighbor;
                parents(child) = parent;
                
                % add child to the list of leaves, so that the list is
                % still sorted
                if isempty(leaves) == 0
                    if degree(child) < degree(leaves(1))
                        resort = 1;
                        index = 0;
                        while resort == 1
                            index = index + 1;
                            if degree(child) > degree(leaves(index))
                                resort = 0;
                                leaves = [leaves(1:index - 1),child,leaves(index:end)];
                            end
                            if index == length(leaves)
                                resort = 0;
                                leaves = [leaves,child];
                            end
                        end
                        
                    else
                        leaves = [child,leaves]; %if degree of child larger than degree of all other leaves
                    end
                else
                    leaves = child; % if no leaves just add to list
                end
                
                % record as found
                tree_size = tree_size + 1;
                node_index(child) = tree_size;
                
                % add edge
                edge_index = edge_index + 1;
                endpoints(edge_index,1) = parent;
                endpoints(edge_index,2) = child;
                
                
                % add to tree
                tree(tree_size - 1) = edge_index;
                parent_edges(child) = edge_index;
                
                
            else
                %% this edge is a chord
                if isnan(searched_from(neighbor)) % make sure we never add a chord twice
                    % add edge
                    edge_index = edge_index + 1;
                    endpoints(edge_index,1) = parent;
                    endpoints(edge_index,2) = child;
                    
                    % add to chords
                    chords_found = chords_found + 1;
                    chords(chords_found) = edge_index;
                    
                end
            end
        end
        
    end
   
    
    %% stop search if no leaves
    if isempty(leaves)
        stop = 1;
        if tree_size < V % check if we've found all nodes
            connected = 0;
        else
            connected = 1;
        end
    end
    
end

%% output
% basic topology
output.node_index = node_index;
output.endpoints = endpoints;

% the tree
output.tree = tree(isnan(tree) == 0);
output.parent_edges = parent_edges;
output.parents = parents;

% the chords
output.chords = chords;

% if connected
output.connected = connected;


end