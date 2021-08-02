function b = getz(A)
% Get MAP class
% Class is defined as {0, 1, 2, ..., (k-1)}
[n, ~] = size(A);
b = zeros(n,1);
for i = 1:n
    [~,temp] = max(A(i,:));
    b(i,1) = temp;
end
b = b-1;
end