function [ fullX ] = fillinBlanks( X, C2C, D2D )
    par.scale = 1;
    par.K = 1;
    par.Top = 5;
    
    [m, n] = size(X);
    if(~iscell(X))
        X = num2cell(X);
    end    

    X_mat = cell2mat(X);
    Pivot_mask = ~(isempty(X_mat) | isnan(X_mat) | isinf(X_mat));        
    Filled_mask = Pivot_mask;
    
    
    for k = 1:par.K
        for c = 1:n
		    unknown_cells = find(cellfun(@(x) isempty(x) | isnan(x) | isinf(x), X(:, c)));
		    filled_cells = setdiff(1:m, unknown_cells);

            %unknown_cells = find(~Pivot_mask(:, c));
            %filled_cells = find(Filled_mask(:, c));
            pivot_cells = find(ismember(filled_cells, find(Pivot_mask(:, c))));
            if(isempty(filled_cells))
                continue;
            end

            for r = 1:numel(unknown_cells)
                row = unknown_cells(r);
                weights = D2D(row, filled_cells);
                if(sum(weights) == 0)
                    continue;
                end
                [~, perm] = sort(weights, 'descend');
                weights(pivot_cells) = par.scale*weights(pivot_cells);
                weights = weights ./ sum(weights);
                
                if(par.Top ~= -1)
                    slim_weight = weights(perm(1:min(numel(weights), par.Top)));
                    slim_filled_cells = filled_cells(perm(1:min(numel(weights), par.Top)));
                else
                    slim_weight = weights;
                    slim_filled_cells = filled_cells;                    
                end
                
                w = 0;
                estim = zeros(size(X(slim_filled_cells(1), c)));
                for i = 1:numel(slim_filled_cells)
                    estim = estim + slim_weight(i)*X{slim_filled_cells(i), c};
		            w = w +  D2D(row, slim_filled_cells(i));
                end
                %estim = estim ./ sum(slim_weight);
%                 estim = estim ./ w;

                X{row, c} = estim;
                Filled_mask(row, c) = true;
            end        
        end



        for r = 1:m
		    unknown_cells = find(cellfun(@(x) isempty(x) | isnan(x) | isinf(x), X(r, :)));
		    filled_cells = setdiff(1:n, unknown_cells);

            %unknown_cells = find(~Pivot_mask(r, :));
            %filled_cells = find(Filled_mask(r, :));

            pivot_cells = find(ismember(filled_cells, find(Pivot_mask(r, :))));

            if(isempty(filled_cells))
                continue;
            end


            for c = 1:numel(unknown_cells)
                col = unknown_cells(c);
                weights = C2C(filled_cells, col);
                if(sum(weights) == 0)
                    continue;
                end            
                [~, perm] = sort(weights, 'descend');
                weights(pivot_cells) = par.scale*weights(pivot_cells);
                weights = weights ./ sum(weights);

                if(par.Top ~= -1)
                    slim_weight = weights(perm(1:min(numel(weights), par.Top)));
                    slim_filled_cells = filled_cells(perm(1:min(numel(weights), par.Top)));
                else
                    slim_weight = weights;
                    slim_filled_cells = filled_cells;                    
                end

                w = 0;
                estim = zeros(size(X(r, slim_filled_cells(1))));
                for i = 1:numel(slim_filled_cells)
                    estim = estim + slim_weight(i)*X{r, slim_filled_cells(i)};
    		        w = w +  C2C(slim_filled_cells(i), col);                    
                end
                % estim = estim ./ sum(slim_weight);
%                 estim = estim ./ w;

                X{r, col} = estim;
                Filled_mask(r, col) = true;            
            end
        end    
    end
    
    fullX = X;
end

