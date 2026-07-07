function result = dot_multiply(A,B)

%求内积<A,B>=tr(AT*B), A,B为矩阵

m = size(A,1);

result = 0;
for i = 1:m
    result = result + A(i,:)*(B(i,:).');
end


end




