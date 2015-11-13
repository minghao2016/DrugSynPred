% runGeneMANIA.m   --> this is a sample script for creating two networks
% and ruuning GeneMANIA using these networks and a list of positive genes. 

%%% the input files are
% Zhang_expData.txt which is microarray data 
% sage_expData_sumTagCounts.txt
%%%%%
%%% 1. import the two files: 

zhang = importdata('Zhang_expData.txt');
sage = importdata('sage_expData_sumTagCounts.txt');
% in these two data files genes are the rows and features are the columns.

%% 2. make two association networks in a cell structure with fields data,
%% collabels, rowlabels --> this is done so that we can use the
%% combineKernel.m function to combine the networks (to make sure they have
%% the same gene order)

K = 50; % set the number of neigbours 

network{1}.data = makeAssociationKernel(zhang.data', K);
network{1}.rowlabels = zhang.textdata(2:end,1);
network{1}.collabels = network{1}.rowlabels;

network{2}.data = makeAssociationKernel(sage.data',K);
network{2}.rowlabels = sage.textdata(2:end,1);
network{2}.collabels = network{2}.rowlabels;

numNetworks = length(network);

%% 3. normalize the networks
for ii = 1:numNetworks
    network{ii}.data = normalizeKernel(network{ii}.data);
end


%% 4. combine the networks
% first create a reference gene list 
geneList = unique([network{1}.rowlabels ; network{2}.rowlabels]);
newNetworks = combineKernels(geneList,network);


%% 5. create a label vector where positive genes are +1, negative genes are
%% -1, and unknown genes are 0. 

%% this is a hypothetical list so we get low area under the AUC curve
[N,N] = size(newNetworks{1}.data); 
labels = zeros(N,1);
posIndx = [1:25];
labels(posIndx) = 1;
labels(labels ~= 1 ) = -1; 

%% 6. get all data field from the structure newNetwork
for ii = 1:numNetworks
    kernels{ii} = newNetworks{ii}.data;
end
    

%% 7. get area under the AUC for prediction, as well as the scores
pp = randperm(N); % this is the permutation index for doing the cross-validation;
nFolds = 3; % number of cross-validation folds
[areas, r, b, weights] = predictClassesCV(labels, kernels, nFolds, pp);


