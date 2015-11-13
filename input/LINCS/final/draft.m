for i = 1:numel(Drugs.Drugs)
max_doses(i) = max(cell2mat(ds3.cdesc(strcmp(ds3.cdesc(:, 7), Drugs.Drugs{i}), 4)));
end


CL_mask = ismember(ds.cdesc(:, 1), CLs);                            
Drug_mask = ismember(ds.cdesc(:, 7), Drugs.Drugs); 
dose_mask = cell2mat(ds.cdesc(:, 4)) == 10;
day_mask = cell2mat(ds.cdesc(:, 9)) == 24; 

final_mask = CL_mask & Type_mask & Drug_mask & dose_mask & day_mask;

ds4 = ds;
ds4.cdesc = ds.cdesc(final_mask, :);  
ds4.cid = ds.cid(final_mask);
ds4.mat = ds.mat(:, final_mask);

CDESC = cell(0, 11);
for i = 1:numel(CLs)
	fprintf('Cellline %s\n', CLs{i});
	for j = 1:numel(Drugs.Drugs)
		fprintf('\tDrug %s\n', Drugs.Drugs{i});
		mask = strcmp(ds4.cdesc(:, 1), CLs{i}) & strcmp(ds4.cdesc(:, 7), Drugs.Drugs{j});
		row = find(mask, 1, 'first');
		fprintf('\t\tNumber of rows = %d\n', numel(row));
		CDESC = [CDESC; ds4.cdesc(row, :)];
	end
end
