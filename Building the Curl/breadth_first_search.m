function output = breadth_first_search(adjacency)
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
tree_size = 1;
chords_found = 0;
current_generation = degree_order(1); % list of leaves of the tree in current generation
node_index(degree_order(1)) = tree_size; 
searched_from = nan(V,1); % equals nan if we haven't searched from this node, 1 otherwise
children = []; % list of children in next generation (new neighbors of current generation)

%% search
stop = 0;
while stop == 0
    
    %% loop over current generation
    for j = 1:length(current_generation)
        parent = current_generation(j);
        neighbors = adjacency{parent};
        searched_from(parent) = 1;
        
        %% loop over neighbors, add to children if not found yet
        for k = 1:length(neighbors)
            neighbor = neighbors(k);
            
            if neighbor ~= parents(parent) % don't backtrack
                if isnan(node_index(neighbor)) % check if already found
                    %% a new node has been found
                    % add to list of children and store parent of child
                    child = neighbor;
                    parents(child) = parent;
                    children = [children,child];
                    
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
                    if isnan(searched_from(neighbor)) % make sure we never add a chord twice, if we'd searched from the neighbor we would have found this edge
                        % add edge
                        edge_index = edge_index + 1;
                        endpoints(edge_index,1) = parent;
                        endpoints(edge_index,2) = neighbor;
                        
                        % add to chords
                        chords_found = chords_found + 1;
                        chords(chords_found) = edge_index;
                        
                    end
                end
            end
        end
    end
    
    %% update current generation
    current_generation = children;
    children = [];
    
    %% cstop search if no new nodes
    if isempty(current_generation)
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