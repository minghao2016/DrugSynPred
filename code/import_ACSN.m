function [ ACSN ] = import_ACSN(  )

% Read vertex annotations
    fd = fopen('input/networks/ACSN/acsn_master_curated.gmt', 'r');
    i = 1;
    while(~feof(fd)) 
        tline = strtrim(fgets(fd));
        tokens = strsplit(tline, '\t');
        ACSN.vertex_names{i, 1} = tokens{1};
        ACSN.vertex_genes{i, 1} = tokens(3:end);
        i = i + 1;
    end
    fclose(fd);
    n = numel(ACSN.vertex_names);

% Read interactome
    EdgeTypes = {'CATALYSIS';'INHIBITION';'MODULATION';'PHYSICAL_STIMULATION';'TRIGGER';'UNKNOWN_CATALYSIS';'UNKNOWN_INHIBITION';'UNKNOWN_POSITIVE_INFLUENCE';'activates';'inhibits';'null'};
    Directed_edges = setdiff(EdgeTypes, 'null');

    fd = fopen('input/networks/ACSN/acsn_ppi.sif', 'r');
    C = textscan(fd, '%[^\t]\t%[^\t]\t%[^\n]\n');
    fclose(fd);

    src = C{1};
    dst = cellfun(@(x) strtrim(x), C{3}, 'UniformOutput', false);
    [~, src_idx] = ismember(src, ACSN.vertex_names);
    [~, dst_idx] = ismember(dst, ACSN.vertex_names);
    directionality_mask = ismember(C{2}, Directed_edges);

    filter_mask = src_idx == 0 | dst_idx == 0;
    src_idx(filter_mask) = [];
    dst_idx(filter_mask) = [];
    directionality_mask(filter_mask) = [];

    ACSN.A = sparse(src_idx, dst_idx, 1, n, n) + sparse(dst_idx(~directionality_mask), src_idx(~directionality_mask), 1, n, n);

% Read functional classes
    fd = fopen('input/networks/ACSN/acsn_v1.1.gmt', 'r');
    i = 1;
    while(~feof(fd)) 
        tline = strtrim(fgets(fd));
        tokens = strsplit(tline, '\t');
        ACSN.class_names{i, 1} = tokens{1};
        ACSN.class_genes{i, 1} = tokens(3:end);
        i = i + 1;
    end
    fclose(fd);
end
