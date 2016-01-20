%% Neural Networks

[~, net] = NNfun(Xtrain, Ytrain_val);

% Make a prediction for the test set
Y_nn = net(Xtest');
Y_nn = round(Y_nn');

% Compute the confusion matrix
C_nn = confusionmat(Ytest_val,Y_nn);
% Examine the confusion matrix for each class as a percentage of the true class
C_nn = bsxfun(@rdivide,C_nn,sum(C_nn,2)) * 100 %#ok<*NOPTS>
sum(diag(C_nn))
%% Generalized Linear Model - Logistic Regression

% Train the classifier
glm = GeneralizedLinearModel.fit(Xtrain,Ytrain_val,'linear','Distribution','binomial','link','logit');

% Make a prediction for the test set
Y_glm = glm.predict(Xtest);
Y_glm = round(Y_glm);

% Compute the confusion matrix
C_glm = confusionmat(Ytest_val,Y_glm);
% Examine the confusion matrix for each class as a percentage of the true class
C_glm = bsxfun(@rdivide,C_glm,sum(C_glm,2)) * 100
sum(diag(C_glm))

%% Discriminant Analysis
da = ClassificationDiscriminant.fit(Xtrain,Ytrain,'discrimType','quadratic');

% Make a prediction for the test set
[Y_da, classifScore_da] = da.predict(Xtest);

% Compute the confusion matrix
C_da = confusionmat(Ytest,Y_da);
% Examine the confusion matrix for each class as a percentage of the true class
C_da = bsxfun(@rdivide,C_da,sum(C_da,2)) * 100
sum(diag(C_da))
%% Classification Using Nearest Neighbors
% Train the classifier
knn = ClassificationKNN.fit(Xtrain,Ytrain,'Distance','seuclidean');

% Make a prediction for the test set
[Y_knn, classifScore_knn] = knn.predict(Xtest);

% Compute the confusion matrix
C_knn = confusionmat(Ytest,Y_knn);
% Examine the confusion matrix for each class as a percentage of the true class
C_knn = bsxfun(@rdivide,C_knn,sum(C_knn,2)) * 100
sum(diag(C_knn))
%% Naive Bayes Classification


% Train the classifier
Nb = NaiveBayes.fit(Xtrain,Ytrain);

% Make a prediction for the test set
Y_Nb = Nb.predict(Xtest);

% Compute the confusion matrix
C_nb = confusionmat(Ytest,Y_Nb);
% Examine the confusion matrix for each class as a percentage of the true class
C_nb = bsxfun(@rdivide,C_nb,sum(C_nb,2)) * 100
sum(diag(C_nb))
%% Support Vector Machines

% Xtrain_temp = Xtrain;
% Xtrain_temp(isnan(Xtrain)) = 0;
% svm_model = svmtrain(2*Ytrain_val-1, zscore(Xtrain_temp), '-s 0 -t 2 -h 0');
% 
% Xtest_temp = Xtest;
% Xtest_temp(isnan(Xtest)) = 0;
% [predicted_label, accuracy, decision_values] = svmpredict(2*Ytest_val-1, zscore(Xtest_temp), svm_model);
% 

opts = statset('MaxIter',30000);
% Train the classifier
svmStruct = svmtrain(Xtrain,Ytrain,'kernel_function','rbf','kktviolationlevel',0.1,'options',opts);

% Make a prediction for the test set
[Y_svm, classifScore_svm]  = svmclassify(svmStruct,Xtest);
C_svm = confusionmat(Ytest,Y_svm);
% Examine the confusion matrix for each class as a percentage of the true class
C_svm = bsxfun(@rdivide,C_svm,sum(C_svm,2)) * 100
sum(diag(C_svm))

%% Decision Trees

tic
% Train the classifier
t = ClassificationTree.fit(Xtrain,Ytrain);
toc

% Make a prediction for the test set
[Y_t, classifScore_ct] = t.predict(Xtest);

% Compute the confusion matrix
C_t = confusionmat(Ytest,Y_t);
% Examine the confusion matrix for each class as a percentage of the true class
C_t = bsxfun(@rdivide,C_t,sum(C_t,2)) * 100
sum(diag(C_t))

%% Ensemble Learning: TreeBagger

% Cost of misclassification
opts = statset('UseParallel',true);
% Train the classifier
tb = TreeBagger(150,Xtrain,Ytrain_val,'method','classification','Options',opts,'OOBVarImp','on','cost',cost);

% Make a prediction for the test set
[Y_tb, classifScore_tb] = tb.predict(Xtest);
Y_tb = nominal(Y_tb);

% Compute the confusion matrix
C_tb = confusionmat(Ytest_val,double(Y_tb)-1);
% Examine the confusion matrix for each class as a percentage of the true class
C_tb = bsxfun(@rdivide,C_tb,sum(C_tb,2)) * 100
sum(diag(C_tb))

%% Compare Results
% This visualization function is making use of a couple files downloaded
% from <http://www.mathworks.com/matlabcentral/ MATLAB Central>, the user
% community website. We are leveraging social computing along the way to
% help us in our effort.

Cmat = [C_nn C_glm C_da C_knn C_nb C_svm C_t C_tb];
labels = {'Neural Net ', 'Logistic Regression ', 'Discriminant Analysis ',...
    'k-nearest Neighbors ', 'Naive Bayes ', 'Support VM ', 'Decision Trees ', 'TreeBagger '};

comparisonPlot( Cmat, labels )
