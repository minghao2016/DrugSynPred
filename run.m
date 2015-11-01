clear
addpath(genpath('code'));

annotations.cellLines = readtable('input/molecular/cell_info.csv', 'Delimiter', ',');
annotations.drugs = readtable('input/synergy/Drugs.csv', 'Delimiter', ',');
annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);


%%
min_threshold = 100;
tissue_names = unique(annotations.cellLines.Tissue__General_);
annotations.groups = {};
for i = 1:numel(tissue_names)
    mask = strcmp( annotations.cellLines.Tissue__General_, tissue_names{i});
    if(nnz(mask) <= min_threshold)
        annotations.groups = [annotations.groups; {annotations.cellLines.Sanger_Name(mask)}];
    else
    end
end

%% Homologous drug identification
tic; WinDTome = readtable('input/DrugHomology/WinDTome.txt', 'Delimiter', '\t', 'Format', '%s %s %d %s %s %s %s %s %s %s %d %d'); toc
[WinDTome_drugs, ic, ia] = unique(WinDTome.Drug_ID);

WinDTome_targets = arrayfun(@(drug_id) unique(WinDTome.Target_Gene_Symbol(ia == drug_id)), 1:numel(WinDTome_drugs), 'UniformOutput', false);

for i = 1:numel(annotations.drugs)
    targets = annotations.drugs.Target{i};
    
%     target_overlap = cellfun(@(D) numel(intersect(D, targets)) / numel(union(D, targets)), WinDTome_targets);
    target_overlap = cellfun(@(D) numel(intersect(D, targets)), WinDTome_targets);

    [~, drug_row, drug_overlap] = find(target_overlap);
    homologous_drug_IDs = WinDTome_drugs(drug_row);
    [~, perm] = sort(zscore(drug_overlap), 'descend');
    homologous_drug_IDs = homologous_drug_IDs(perm);
    
end

% Now match IC50
X = readtable('input/DrugHomology/gdsc_manova_input_w5.csv');
