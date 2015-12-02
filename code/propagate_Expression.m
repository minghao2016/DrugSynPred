
function [D2D_prop_Exp,C2C_prop_Exp] = propagate_Expression(D2D, C2C , Expr_DS, alpha)

    
    [numberofDrug numberofCell] = size(Expr_DS);
    %number of gene expression
    numberofGenes = length(Expr_DS{find(~cellfun(@isempty,Expr_DS),1)});
    
    fprintf('Propagating expression through drug similarity network ...\n');
    %Propagate through D2D
    geneExpression = zeros(numberofGenes,numberofDrug);
    for i=1:numberofCell
        known = find(~cellfun(@isempty,Expr_DS(:,i)));
        if(known)
            geneExpression(:,known) = cell2mat(Expr_DS(known,i)');
            D2D_prop_Exp = propagation(D2D,alpha,geneExpression');
        end
    end
    
    fprintf('Propagating expression through cell similarity network ...\n');
    %Propagate through C2C
    geneExpression = zeros(numberofGenes,numberofCell);
    for i=1:numberofDrug
        known = find(~cellfun(@isempty,Expr_DS(i,:)));
        if(known)
            geneExpression(:,known) = cell2mat(Expr_DS(i,known));
            C2C_prop_Exp = propagation(C2C,alpha,geneExpression');
        end
    end
    
end
    
    
    