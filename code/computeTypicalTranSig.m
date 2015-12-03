function [ func_synergy_weight, CC, x_report] = computeTypicalTranSig( varargin  )
    params = inputParser;
    params.addParamValue('visualize', false, @(x) islogical(x) );
    params.parse(varargin{:});
    par = params.Results;
    
    load('/home/shahin/Dropbox/Dream/experiment/input/LINCS/final/LINCS_Class_gene_expressions_Z2.mat');
    x = cell2mat(cellfun(@(class_expr) quantile(sum(class_expr, 1), [.05, 0.5, 0.95]) ./ std(sum(class_expr, 1)), XXX, 'UniformOutput', false));
    delta = x(:, 3) - abs(x(:, 1));
    func_synergy_weight = delta ./ norm(delta, 1);
    x_report = [ACSN.class_names, num2cell(100*func_synergy_weight)];
    [~, perm] = sort(delta);
    x_report = x_report(perm, :);
    
    X = cell2mat(cellfun(@(class_expr) sum(class_expr, 1) ./ std(sum(class_expr, 1)), XXX, 'UniformOutput', false))';
    [CC] = corr(X);
    if(par.visualize)
        clustergram(-CC, 'RowLabels', ACSN.class_names, 'ColumnLabels', ACSN.class_names, 'Linkage', 'average', 'OPTIMALLEAFORDER', true, 'Standardize', 1);
    end
end

