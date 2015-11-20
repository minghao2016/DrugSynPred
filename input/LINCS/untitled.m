CL_mask = ismember(ds.cdesc(:, 1), CLs);                            
Drug_mask = ismember(ds.cdesc(:, 7), Drugs.Drugs); 
dose_mask = cell2mat(ds.cdesc(:, 4)) == 10;
day_mask = cell2mat(ds.cdesc(:, 9)) == 24; 

final_mask = CL_mask & Type_mask & Drug_mask & dose_mask & day_mask;

ds_sub = ds;
ds_sub.cdesc = ds.cdesc(final_mask, :);  
ds_sub.cid = ds.cid(final_mask);
ds_sub.mat = ds.mat(:, final_mask);

% 152/180 (36*5) of pairs
idx = 1;
CDESC = {};
for i = 1:numel(CLs)
	fprintf('Cellline %s\n', CLs{i});
	for j = 1:numel(Drugs.Drugs)
		fprintf('\tDrug %s\n', Drugs.Drugs{i});
		mask = strcmp(ds_sub.cdesc(:, 1), CLs{i}) & strcmp(ds_sub.cdesc(:, 7), Drugs.Drugs{j});
		rows = find(mask);
        if(numel(rows) == 0)
            continue;
        end
		fprintf('\t\tNumber of rows = %d\n', numel(rows));
		first_row = find(mask, 1, 'first');        
		CDESC(idx, 1:11) = ds_sub.cdesc(first_row, :);
        MAT(1:22268, idx) = mean(ds_sub.mat(:, rows), 2);
        CID(idx) = ds_sub.cid(first_row);
        idx = idx + 1;
	end
end

final_ds = ds_sub;
final_ds.cdesc = CDESC;
final_ds.mat = MAT;
final_ds.cid = CID;

mkgct('final_ds.gct', final_ds)
