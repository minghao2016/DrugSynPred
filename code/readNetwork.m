function [ Net ] = readNetwork(  )
        % Expand network using ESEA edges
        fd = fopen('input/networks/ESEA/ESEA_pathway_edges.txt', 'r');
        C = textscan(fd, '%s %s');
        fclose(fd);
        Net.vertex_names = union(C{1}, C{2});
        Net.vertex_genes = Net.vertex_names;        
        
        [~, ii] = ismember(C{1}, Net.vertex_names);
        [~, jj] = ismember(C{2}, Net.vertex_names);
        Net.A = sparse(ii, jj, 1, numel(Net.vertex_names), numel(Net.vertex_names));
        Net.A = max(Net.A, Net.A');
        
        % We can combine this with SignaLink to add directions, especially
        % for signaling edges, as well as GRN


end

