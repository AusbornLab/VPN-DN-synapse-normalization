
close all
clear all



start_trees
load_tree()

% List of synapse data Excel files
file_names = {
    'hemibrain_DNp01_inputs.xlsx',
};

% VPN types of interest
VPNs = {'LC4','LC6', 'LC22', 'LPLC1', 'LPLC2', 'LPLC4'};

% Color map for each VPN
vpn_colors = containers.Map( ...
    {'LC4','LC6','LC22','LPLC1','LPLC2','LPLC4'}, ...
    {[0 0 1], [1, 0, 0.91], [1 1 0], [1 0 0], [1 0.647 0], [0 1 0]} ...
);

% Loop through each file
for i = 1:length(file_names)
    data = readtable(file_names{i});
    
    % Filter by VPN type and sort
    data_filtered = data(ismember(data.type, VPNs), :);
    data_filtered = sortrows(data_filtered, 'type');

    % Extract XYZ coordinates
    syn_coords = [data_filtered.post_x, data_filtered.post_y, data_filtered.post_z];

    % Match each synapse to nearest node in the tree
    node_coords = [trees{1}.X, trees{1}.Y, trees{1}.Z];
    synapseNodes = knnsearch(node_coords, syn_coords);

    % Prepare for plotting dendrogram
    [xdend, tree_dend] = xdend_tree(trees{1});
yvec = Pvec_tree(trees{1});

% Create parent-child edge indices
idpar = trees{1}.dA * (1:size(trees{1}.dA, 1))';
idpar(idpar == 0) = 1;

% Original edge coordinates
X1 = xdend(idpar); 
X2 = xdend;
Y1 = yvec(idpar);  
Y2 = yvec;

% Expanded arrays to explicitly define start and end points of each edge
X1all = [X1; X2];
X2all = [X2; X2];
Y1all = [Y1; Y1];
Y2all = [Y1; Y2];

% Plot dendrogram lines
figure;
line([Y1all Y2all]', [X1all X2all]', 'Color', 'k'); 
hold on;

% Plot each VPN type with color and synapse markers
for vpn = VPNs
    vpn_name = vpn{1};
    rgb = vpn_colors(vpn_name);
    indices = strcmp(data_filtered.type, vpn_name);
    nodes = synapseNodes(indices);
    plot(yvec(nodes), xdend(nodes), 'o', ...
        'Color', rgb, 'MarkerFaceColor', rgb, 'MarkerSize', 3);
end

title(sprintf('Dendrogram with Synapses: %s', file_names{i}), 'Interpreter', 'none');
xlabel('Dendritic path length');
ylabel('Branch order');
view(0, 90);
axis xy;
    ylim([-20 800])
    xlim([0 350])
set(gcf, 'Position', [100 100 900 500]);
end

