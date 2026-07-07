function [U, V, iter, time] = V_cycle(F, G, v1, v2, n, l, MAXITER)
%用V-cycle方法解原Stokes方程


%v1,v2:下降上升过程中分别的迭代次数
%n:log2问题规模
%l:log2底层维度
%MAXITER:最大V-cycle迭代次数
%iter:进行的V-cycle次数

tic;

N = 2^n;
iter1 = n-l;

r0 = sqrt(norm(F,'fro')^2 + norm(G,'fro')^2);

%存储每一层的U,V,P,F,G,Z
U_store = cell(1,iter1+1);
V_store = cell(1,iter1+1);
P_store = cell(1,iter1+1);
F_store = cell(1,iter1);
G_store = cell(1,iter1);
Z_store = cell(1,iter1);

U_store{1} = zeros(N+1,N);
V_store{1} = zeros(N,N+1);
P_store{1} = zeros(N,N);
F_store{1} = F;
G_store{1} = G;
Z_store{1} = zeros(N,N);


iter = 0;
rh = r0;

while rh/r0 > 1e-8 && iter < MAXITER
    iter = iter+1;
    %下降过程
    for i = 1:iter1
        if i ~= 1
            %初值为0
            U = zeros(N+1,N);
            V = zeros(N,N+1);
            P = zeros(N,N);
        else
            U = U_store{1};
            V = V_store{1};
            P = P_store{1};
            F = F_store{1};
            G = G_store{1};
            Z = Z_store{1};
        end

        %用v1次DGS迭代求解得X(v1)
        for k = 1:v1
            [U,V,P] = DGS(U,V,P,F,G,Z,N);
        end

        %记录这一层的U,V,P,F,G,Z
        U_store{i} = U;
        V_store{i} = V;
        P_store{i} = P;

        if i ~= 1 
            F_store{i} = F;
            G_store{i} = G;
            Z_store{i} = Z;
        end
        
        %计算残差r, 应用限制算子进入下一层
        [U1,V1,P1] = left_multiply_matrix(U,V,P,N);
        [F,G,Z] = limit(F-U1,G-V1,Z-P1,N);
        N = N/2;
         
    end
        

    %此时到达最底层, 将DGS迭代解近似作为真解

    U = zeros(N+1,N);
    V = zeros(N,N+1);
    P = zeros(N,N);

    for k = 1:100
        [U,V,P] = DGS(U,V,P,F,G,Z,N);
    end

    U_store{iter1+1} = U;
    V_store{iter1+1} = V;
    P_store{iter1+1} = P;

    %上升过程

    for i = iter1:-1:1
        %对r用提升算子, 从i+1层提到i层
        [U,V,P] = lift(U,V,P,N);
        N = N*2;

        %计算第i层的X(v1+1)
        U = U_store{i} + U;
        V = V_store{i} + V;
        P = P_store{i} + P;

        %再进行v2次DGS迭代,得到第i层的新解X(v1+v2+1)
        for k = 1:v2
            [U,V,P] = DGS(U,V,P,F_store{i},G_store{i},Z_store{i},N);
        end

        if i == 1
            U_store{1} = U;
            V_store{1} = V;
            P_store{1} = P;
        end
        
    end

    %计算rh/r0
    [U1,V1,P1] = left_multiply_matrix(U,V,P,N);
    rh = sqrt(norm(U1-F_store{1})^2 + norm(V1-G_store{1})^2 + norm(P1,2)^2);

    disp(['V-cycle',num2str(iter),': rh=',num2str(rh)]);

end

time = toc;


end
