function P1 = left_multiply_BT(U, V, N)

%实际上是把矩阵U,V写成向量形式后,左乘BT


h = 1/N;


P1 = zeros(N,N);

for i = 1:N
    for j = 1:N
        P1(i,j) = -1/h*(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j));
    end
end


end