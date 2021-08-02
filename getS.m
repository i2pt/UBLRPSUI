function sum = getS(A,b)
% size(A) = [n x p]
% size(b) = [p x 1]
[n, p] = size(A);
if length(b) ~= p
    disp('Error : Check Matrix & Vector Size <mahayat>')
end
sum = 0;
for i = 1:n
    sum = sum + (A(i,:)'-b)*((A(i,:)'-b))';
end