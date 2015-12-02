addpath(genpath('lib'))
ds = parse_gctx('GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n78980x22268_2015-06-30.gct');

CLs = {'A549', 'BT20', 'HS578T', 'HT29', 'MCF7', 'MDAMB231'};
Drugs = {'AZD-7762', 'AZD-8055', 'AZD4547', 'BYL719', 'CYT387', 'GDC-0941', 'GSK-1059615', 'HG-14-10-04', 'KU-55933', 'MK-1775', 'MK-2206', 'motesanib', 'NU-7441', 'NVP-AUY922', 'NVP-TAE684', 'PD-173074', 'PF-04217903', 'azacitidine', 'celastrol', 'crizotinib', 'entinostat', 'enzastaurin', 'everolimus', 'gefitinib', 'lapatinib', 'linsitinib', 'neratinib', 'olaparib', 'regorafenib', 'saracatinib', 'selumetinib', 'sorafenib', 'tivantinib', 'trametinib', 'vandetanib', 'vorinostat'};
all_Drugs = {'(+)-JQ1', '(-)-JQ1', 'A 769662', 'A443654', 'A66', 'ABT-737', 'AGK-2', 'AKT-inhibitor-1-2', 'ALW-II-38-3', 'ALW-II-49-7', 'AS-601245', 'AS-605240', 'AST1306', 'AT-7519', 'AT7867', 'AZ-628', 'AZ20', 'AZD-1480', 'AZD-5438', 'AZD-6482', 'AZD-7762', 'AZD-8055', 'AZD-8330', 'AZD4547', 'AZD7762', 'BI-2536', 'BIX-01294', 'BIX-02189', 'BMS-345541', 'BMS-387032', 'BMS-777607', 'BMS509744', 'BS-181', 'BX-795', 'BX-912', 'BYL719', 'C646', 'CC-401', 'CG-930', 'CGK-733', 'CGP-60474', 'CH5424802', 'CHIR-99021', 'CP-724714', 'CP466722', 'CTB', 'CX-5461', 'CYT387', 'D-4476', 'DCC-2036', 'DMSO', 'EX-527', 'FR-180204', 'GDC-0068', 'GDC-0879', 'GDC-0941', 'GDC-0980', 'GSK-1059615', 'GSK-1070916', 'GSK-1904529A', 'GSK-2126458', 'GSK-2334470', 'GSK-3-inhibitor-II', 'GSK-429286A', 'GSK-461364', 'GSK-690693', 'GSK-J1', 'GSK-J2', 'GSK-J4', 'GSK1059615', 'GW-5074', 'GW-843682X', 'HG-14-10-04', 'HG-14-8-02', 'HG-5-113-01', 'HG-5-88-01', 'HG-6-64-01', 'HG-9-91-01', 'HY-11007', 'I-BET', 'I-BET151', 'INK-128', 'IOX2', 'JNJ-38877605', 'JNK-9L', 'JNK-IN-5A', 'JW-7-24-1', 'JW55', 'JWE-035', 'KI20227', 'KIN001-043', 'KIN001-055', 'KIN001-127', 'KIN001-220', 'KIN001-242', 'KIN001-244', 'KIN001-260;', 'KIN001-265', 'KIN001-266', 'KIN001-269', 'KIN001-270', 'KIN236', 'KU-55933', 'KU-60019', 'KW-2449', 'LDN-193189', 'LY2603618', 'MC1568', 'MGCD-265', 'MK-1775', 'MK-2206', 'MLN-8054', 'NPK76-II-72-1', 'NU-7026', 'NU-7441', 'NVP-AEW541', 'NVP-AUY922', 'NVP-BEZ235', 'NVP-BGJ398', 'NVP-BGT226', 'NVP-TAE226', 'NVP-TAE684', 'ON-01910', 'OSI-027', 'OSI-930', 'OTSSP167', 'PD-0325901', 'PD-173074', 'PD-184352', 'PF 573228', 'PF-04217903', 'PF-3758309', 'PF-431396', 'PF-477736', 'PF-562271', 'PFI-1', 'PHA-665752', 'PHA-767491', 'PHA-793887', 'PI-103', 'PIK-93', 'PLX-4720', 'QL-X-138', 'QL-XI-92', 'QL-XII-47', 'RAF 265', 'RG-108', 'SAR-245408', 'SAR-245409', 'SB-203580', 'SB-216763', 'SB-239063', 'SB-525334', 'SB590885', 'SGI-1776', 'SGX523', 'SU-11274', 'SYK-inhibitor', 'TAK-715', 'TG-101348', 'TGX-221', 'THZ-2-98-01', 'TPCA-1', 'TWS-119', 'UNC0638', 'UNC1215', 'UNC669', 'VE821', 'VX-745', 'WH-4-025', 'WYE-125132', 'WZ-3105', 'WZ-4-145', 'WZ-4002', 'WZ-7043', 'XMD-1150', 'XMD-1499', 'XMD-885', 'XMD-892', 'XMD11-85H', 'XMD16-144', 'Y-27632', 'Y-39983', 'YM-201636', 'ZM-447439', 'ZSTK-474', 'afatinib', 'alvocidib', 'amuvatinib', 'anacardic-acid', 'azacitidine', 'barasertib', 'belinostat', 'bosutinib', 'brivanib', 'buparlisib', 'cabozantinib', 'canertinib', 'celastrol', 'chelerythrine chloride', 'crizotinib', 'dabrafenib', 'dacomitinib', 'daminozide', 'dasatinib', 'decitabine', 'dinaciclib', 'doramapimod', 'dovitinib', 'entinostat', 'enzastaurin', 'epigallocatechin', 'erlotinib', 'everolimus', 'foretinib', 'fostamatinib', 'garcinol', 'gefitinib', 'geldanamycin', 'ibrutinib', 'idelalisib', 'imatinib', 'iniparib', 'ischemin', 'lapatinib', 'linifanib', 'linsitinib', 'masitinib', 'methylstat', 'mitoxantrone', 'mocetinostat', 'motesanib', 'neratinib', 'nilotinib', 'nintedanib', 'olaparib', 'palbociclib', 'pazopanib', 'pelitinib', 'ponatinib', 'pracinostat', 'pyrazolanthrone', 'quizartinib', 'radicicol', 'regorafenib', 'resveratrol', 'rocilinostat', 'rondual-kinase-inhibitor', 'roscovitine', 'rucaparib', 'ruxolitinib', 'saracatinib', 'selumetinib', 'serdemetan', 'sirolimus', 'sorafenib', 'sunitinib', 'sutent', 'tanespimycin', 'tivantinib', 'tivozanib', 'torin-1', 'torin-2', 'tozasertib', 'trametinib', 'tranylcypromine', 'tyrphostin-AG-1478', 'vandetanib', 'veliparib', 'vemurafenib', 'vorinostat', 'withaferin-a'};


CL_mask = ismember(ds.cdesc(:, 1), CLs);                            
Drug_mask = ismember(ds.cdesc(:, 7), Drugs); 

dose_mask = cell2mat(ds.cdesc(:, 4)) == 10;
day_mask = cell2mat(ds.cdesc(:, 9)) == 24; 

final_mask = CL_mask & Drug_mask & dose_mask & day_mask;

ds_sub = ds;
ds_sub.cdesc = ds.cdesc(final_mask, :);  
ds_sub.cid = ds.cid(final_mask);
ds_sub.mat = ds.mat(:, final_mask);

% 152/180 (36*5) of pairs
idx = 1;
CDESC = {};
for i = 1:numel(CLs)
	fprintf('Cellline %s\n', CLs{i});
	for j = 1:numel(Drugs)
		fprintf('\tDrug %s\n', Drugs{i});
		mask = strcmp(ds_sub.cdesc(:, 1), CLs{i}) & strcmp(ds_sub.cdesc(:, 7), Drugs{j});
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



MAT = ds.mat(:, CL_mask & dose_mask & day_mask);
MAT(MAT > -1.65 & MAT < 1.65) = 0; 
ACSN_class_rows = cellfun(@(genes) find(ismember(ds.rdesc(:, 7), genes)), ACSN.class_genes, 'UniformOutput', false);

