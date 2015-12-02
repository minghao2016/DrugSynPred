
function [Expr_DS_Prop_D2D,Expr_DS_Prop_C2C] = propagate_Expression(D2D, C2C , Expr_DS, alpha)

    
    [numberofDrug numberofCell] = size(Expr_DS);
    %number of gene expression
    numberofGenes = length(Expr_DS{find(~cellfun(@isempty,Expr_DS),1)});
    
    fprintf('Propagating expression through drug similarity network ...\n');
    %Propagate through D2D
    geneExpression = zeros(numberofGenes,numberofDrug);
    Expr_DS_Prop_D2D=cell(numberofDrug, numberofCell);
    for i=1:numberofCell
        known = find(~cellfun(@isempty,Expr_DS(:,i)));
        if(known)
            geneExpression(:,known) = cell2mat(Expr_DS(known,i)');
            D2D_prop_Exp = propagation(D2D,alpha,geneExpression');
            
            for j=1:numberofDrug
                Expr_DS_Prop_D2D{j,i} = num2cell(D2D_prop_Exp(j,:));
            end
        end
       
    end
    
    fprintf('Propagating expression through cell similarity network ...\n');
    %Propagate through C2C
    geneExpression = zeros(numberofGenes,numberofCell);
    Expr_DS_Prop_C2C=cell(numberofDrug, numberofCell);
    for i=1:numberofDrug
        known = find(~cellfun(@isempty,Expr_DS(i,:)));
        if(known)
            geneExpression(:,known) = cell2mat(Expr_DS(i,known));
            C2C_prop_Exp = propagation(C2C,alpha,geneExpression');
            for j=1:numberofCell
                Expr_DS_Prop_C2C{i,j} = num2cell(C2C_prop_Exp(j,:));
            end
        end
         
    end
    
end
    
    
    