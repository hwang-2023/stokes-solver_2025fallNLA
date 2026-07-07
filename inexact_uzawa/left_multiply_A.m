function [U1, V1] = left_multiply_A(U, V, N)

%实际上是把矩阵U,V写成向量形式后,左乘A


h = 1/N;

U1 = zeros(N+1,N);

for j = 2:N-1
    for i = 2:N
        U1(i,j) = -1/(h^2)*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)-4*U(i,j));
    end
end

for i = 2:N
    U1(i,1) = -1/(h^2)*(U(i-1,1)+U(i+1,1)-2*U(i,1)+U(i,2)-U(i,1));
    U1(i,N) = -1/(h^2)*(U(i-1,N)+U(i+1,N)-2*U(i,N)-U(i,N)+U(i,N-1));
end

U1(1,1:N) = U(1,1:N);
U1(N+1,1:N) = U(N+1,1:N);


V1 = zeros(N,N+1);


for i = 2:N-1
    for j = 2:N
        V1(i,j) = -1/(h^2)*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j));
    end
end

for j = 2:N
    V1(1,j) = -1/(h^2)*(V(2,j)-V(1,j)+V(1,j+1)+V(1,j-1)-2*V(1,j));
    V1(N,j) = -1/(h^2)*(V(N-1,j)-V(N,j)+V(N,j+1)+V(N,j-1)-2*V(N,j));
end

V1(1:N,1) = V(1:N,1);
V1(1:N,N+1) = V(1:N,N+1);


end




