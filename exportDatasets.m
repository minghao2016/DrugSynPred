    clear
    load datasets
    D2D = Construct_D2D(annotations, interactome);
    [C2C, NodeWeights] = Construct_C2C(annotations, interactome, 'expression_only', true);
    
   
    X = cell2mat(fillinBlanks(Mono.EMax, C2C, D2D));    
    X = X ./ 100;

    X_ds = mat2dataset(X, 'VarNames', annotations.cellLines.Sanger_Name, 'ObsNames', annotations.drugs.ChallengeName);    
    export(X_ds, 'file', 'Dream_sensitivity.csv', 'Delimiter', ',');
    
    annotations.drugs.Target
    
    all_targets = {};
    for i = 1:size(annotations.drugs, 1)
        all_targets = union(all_targets, annotations.drugs.Target{i});
    end
    
    Y = cell2mat(arrayfun(@(d) ismember(all_targets, annotations.drugs.Target{d})', 1:size(annotations.drugs, 1), 'UniformOutput', false))';
    Y_ds = mat2dataset(Y, 'VarNames', all_targets, 'ObsNames', annotations.drugs.ChallengeName);    
    export(Y_ds, 'file', 'Dream_interaction_binary.csv', 'Delimiter', ',');
    
    
    