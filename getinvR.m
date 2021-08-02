function invR = getinvR(R)
[a, b] = size(R);
m = a/b;
invR = zeros(a,b);
for i = 0:m-1
    t = ((i*b)+1):(b*(i+1));
    invR(t,:) = inv(R(t,:));
end
end