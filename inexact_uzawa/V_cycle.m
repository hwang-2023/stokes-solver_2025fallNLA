function [U, V] = V_cycle(F, G, v1, v2, N, L)

%求解PCG中的MZ=R
%实际是V-cycle法近似求解AZ=R

%v1,v2:下降上升过程中分别的迭代次数
%N:问题规模

%这里取最底层为L=4
iter1 = log2(N)-log2(L);


U_store = cell(1,iter1+1);
V_store = cell(1,iter1+1);

F_store = cell(1,iter1);
G_store = cell(1,iter1);


U_store{1} = zeros(N+1,N);
V_store{1} = zeros(N,N+1);

F_store{1} = F;
G_store{1} = G;

r0 = sqrt(norm(F,'fro')^2 + norm(G,'fro')^2);
rh = r0;
iter = 0;

while rh/r0 > 1e-4 && iter < 5
    iter = iter + 1;
    %下降过程
    for i = 1:iter1
        if i ~= 1
            %初值为0
            U = zeros(N+1,N);
            V = zeros(N,N+1);
        else
            U = U_store{1};
            V = V_store{1};
            F = F_store{1};
            G = G_store{1};
        end

        %用v1次G-S迭代求解得X(v1)
        for k = 1:v1
            [U,V] = GS(U,V,F,G,N);
        end

        %记录这一层的U,V,F,G
        U_store{i} = U;
        V_store{i} = V;

        if i ~= 1 
            F_store{i} = F;
            G_store{i} = G;
        end
        
        %计算r,限制算子进入下一层
        [U1,V1] = left_multiply_A(U,V,N);
        [F,G] = limit(F-U1,G-V1,N);
        N = N/2;
         
    end
            U = zeros(N+1,N);
            V = zeros(N,N+1);

    %此时到达最底层,将G-S迭代解近似作为真解
        for k = 1:100
            [U,V] = GS(U,V,F,G,N);
        end

    U_store{iter1+1} = U;
    V_store{iter1+1} = V;


    %上升过程

    for i = iter1:-1:1
        %对r用提升算子,从i+1层提到i层
        [U,V] = lift(U_store{i+1},V_store{i+1},N);
        N = N*2;

        %计算第i层的X(v1+1)
        U = U_store{i} + U;
        V = V_store{i} + V;

        %再进行v2次G-S迭代,得到第i层的新解X(v1+v2+1)
        for k = 1:v2
            [U,V] = GS(U,V,F_store{i},G_store{i},N);
        end
        U_store{i} = U;
        V_store{i} = V;
    end

    [U1,V1] = left_multiply_A(U,V,N);
    rh = sqrt(norm(U1-F_store{1})^2 + norm(V1-G_store{1})^2);

end

end







