function [ fullX ] = fillinBlanks( X, C2C, D2D )
    [m, n] = size(X);
    if(~iscell(X))
        X = num2cell(X);
    end    
    
    % fill in columns
    for c = 1:n
        empty_cells = find(cellfun(@(x) isempty(x) | isnan(x) | isinf(x), X(:, c)));
        filled_cells = setdiff(1:m, empty_cells);
        if(isempty(filled_cells))
            continue;
        end
            
        weights = D2D(row, filled_cells);
        weights = weights ./ sum(weights);
        for r = 1:numel(empty_cells)
            row = empty_cells(r);
            estim = zeros(size(X(filled_cells(1), c)));
            for i = 1:numel(filled_cells)
                estim = estim + weights(i)*X{filled_cells(i), c};
            end
        end
        estim = estim ./ sum(weights);
    end

    for r = 1:m
        empty_cells = find(cellfun(@(x) isempty(x) | isnan(x) | isinf(x), X(r, :)));
        filled_cells = setdiff(1:n, empty_cells);
        if(isempty(filled_cells))
            continue;
        end
            
        weights = C2C(filled_cells, col);
        weights = weights ./ sum(weights);
        for c = 1:numel(empty_cells)
            col = empty_cells(c);
            estim = zeros(size(X(r, filled_cells(1))));
            for i = 1:numel(filled_cells)
                estim = estim + weights(i)*X{r, filled_cells(i)};
            end
        end
        estim = estim ./ sum(weights);
    end    
    fullX = X;
end

