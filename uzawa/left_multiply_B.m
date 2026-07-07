function [U1, V1] = left_multiply_B(P, N)

%实际上是把矩阵P写成向量形式后,左乘B


h = 1/N;
U1 = zeros(N+1,N);
V1 = zeros(N,N+1);


for j = 2:N-1
    for i = 2:N
        U1(i,j)=  1/h*(P(i,j)-P(i-1,j));
    end
end

for i = 2:N
    U1(i,1) = 1/h*(P(i,1)-P(i-1,1));
    U1(i,N) = 1/h*(P(i,N)-P(i-1,N));
end

for i = 2:N-1
    for j = 2:N
        V1(i,j) = 1/h*(P(i,j)-P(i,j-1));
    end
end

for j = 2:N
    V1(1,j) = 1/h*(P(1,j)-P(1,j-1));
    V1(N,j) = 1/h*(P(N,j)-P(N,j-1));
end


end



