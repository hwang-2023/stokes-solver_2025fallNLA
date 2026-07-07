function [U1, V1, P1] = left_multiply_matrix(U, V, P, N)

%实际上是把矩阵U,V,P写成向量形式后,左乘
% | A  B | 
% | BT 0 |


h = 1/N;

U1 = zeros(N+1,N);

for j = 2:N-1
    for i = 2:N
        U1(i,j) = -1/(h^2)*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)-4*U(i,j));
        U1(i,j) = U1(i,j) + 1/h*(P(i,j)-P(i-1,j));
    end
end

for i = 2:N
    U1(i,1) = -1/(h^2)*(U(i-1,1)+U(i+1,1)-2*U(i,1)+U(i,2)-U(i,1));
    U1(i,1) = U1(i,1) + 1/h*(P(i,1)-P(i-1,1));
    U1(i,N) = -1/(h^2)*(U(i-1,N)+U(i+1,N)-2*U(i,N)-U(i,N)+U(i,N-1));
    U1(i,N) = U1(i,N) + 1/h*(P(i,N)-P(i-1,N));
end

U1(1,1:N) = U(1,1:N);
U1(N+1,1:N) = U(N+1,1:N);


V1 = zeros(N,N+1);


for i = 2:N-1
    for j = 2:N
        V1(i,j) = -1/(h^2)*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j));
        V1(i,j) = V1(i,j) + 1/h*(P(i,j)-P(i,j-1));
    end
end

for j = 2:N
    V1(1,j) = -1/(h^2)*(V(2,j)-V(1,j)+V(1,j+1)+V(1,j-1)-2*V(1,j));
    V1(1,j) = V1(1,j) + 1/h*(P(1,j)-P(1,j-1));
    V1(N,j) = -1/(h^2)*(V(N-1,j)-V(N,j)+V(N,j+1)+V(N,j-1)-2*V(N,j));
    V1(N,j) = V1(N,j) + 1/h*(P(N,j)-P(N,j-1));
end

V1(1:N,1) = V(1:N,1);
V1(1:N,N+1) = V(1:N,N+1);

P1 = zeros(N,N);

for i = 1:N
    for j = 1:N
        P1(i,j) = -1/h*(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j));
    end
end


end




