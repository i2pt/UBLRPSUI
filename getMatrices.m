function [TPR, PRE, F1, ACC] = getMatrices(A)
TN = A(1,1);
FP = A(1,2);
FN = A(2,1);
TP = A(2,2);

TPR = TP/(TP+FN);
TNR = TN/(TN+FP);
PRE = TP/(TP+FP);
F1 = 2/((1/TPR)+(1/PRE));
ACC = (TP+TN)/(TP+TN+FP+FN);
end