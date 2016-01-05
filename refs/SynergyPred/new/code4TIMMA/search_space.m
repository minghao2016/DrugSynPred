function [IM_d,IM_i,IM_o] = search_space(ds, k_set, S, y_vals)
% ds = number of drugs
% k_set = the current kinase set
% S = the complete data

a = sum(k_set)+1;
[rows,cols,G_dec] = graycode2(a);

IM_d = NaN(rows,cols,ds,2);
IM_i = NaN(rows,cols,ds,2);
IM_o = NaN(rows,cols,ds,2);

for i = 1:ds
    tmp1 = [S(i,k_set==1) 0];
   dec_val= bin2dec(char(tmp1 +48));
%     tmp111 = double(cellfun(@(x) all(x==tmp1),G));
     tmp11 = G_dec==dec_val;
     tmp11 = +tmp11;
%    
    [super_set_0,sub_set_0]=binary_set([S(i,k_set==1) 0]);
    tmp12=+ismember(G_dec,super_set_0);
    tmp13=+ismember(G_dec,sub_set_0);
    tmp11(tmp11==0)=NaN;
    tmp12(tmp12==0)=NaN;
    tmp13(tmp13==0)=NaN;
% 
    IM_d(:,:,i,1) = tmp11.*y_vals(i);
    IM_i(:,:,i,1) = tmp12.*y_vals(i);
    IM_o(:,:,i,1) = tmp13.*y_vals(i);
    
    tmp2 = [S(i,k_set==1) 1];
     dec_val2= bin2dec(char(tmp2 +48));
      tmp21 = G_dec==dec_val2;
     tmp21 = +tmp21;
    [super_set_1,sub_set_1]=binary_set(tmp2);
    tmp22=+ismember(G_dec,super_set_1);
    tmp23=+ismember(G_dec,sub_set_1);
%     tmp211 = double(cellfun(@(x) all(x==tmp2),G));
%     tmp222 = double(cellfun(@(x) all(x(tmp2==1)),G));
%     tmp233 = double(~tmp211.*cellfun(@(x) all(tmp2(x==1)),G));
%     
%      tmp211(tmp211==0)=NaN;
%     tmp222(tmp222==0)=NaN;
%     tmp233(tmp233==0)=NaN;
    tmp21(tmp21==0)=NaN;
    tmp22(tmp22==0)=NaN;
    tmp23(tmp23==0)=NaN;

    IM_d(:,:,i,2) = tmp21.*y_vals(i);
    IM_i(:,:,i,2) = tmp22.*y_vals(i);
    IM_o(:,:,i,2) = tmp23.*y_vals(i);
end
