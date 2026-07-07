function [U, V, iter] = PCG(F, G, N, L, v1, v2, epsilon)


    %V-cycle做预条件的预优共轭梯度法求解AU = F - BP
    %F,G:RHS的矩阵形式
    %epsilon:停止条件tau*||BTU||

    b = sqrt(norm(F,'fro')^2+norm(G,'fro')^2);
    rho = b^2;
    delta = 10*epsilon;
    iter = 0;

    U = zeros(N+1,N);
    V = zeros(N,N+1);

    RU = F;
    RV = G;

    %由于预优共轭梯度法收敛较快, 取一个较小的最大迭代次数
    MAXITER = N;

    while sqrt(rho) > b*1e-8 && delta > epsilon && iter < MAXITER
        iter = iter + 1;
        %Mz=r
        [ZU,ZV] = V_cycle(RU,RV,v1,v2,N,L);
        if iter == 1
            PU = ZU;
            PV = ZV;
            rho = dot_multiply(RU,ZU) + dot_multiply(RV,ZV);
        else
            rho_before = rho;
            %rho = <z,r>
            rho = dot_multiply(RU,ZU) + dot_multiply(RV,ZV);
            beta = rho/rho_before;
            %p=z+beta*p
            PU = ZU + beta*PU;
            PV = ZV + beta*PV;
        end
        [WU,WV] = left_multiply_A(PU,PV,N);
        alpha = rho/(dot_multiply(WU,PU)+dot_multiply(WV,PV));
        %u = u+alpha*u
        U = U + alpha*PU;
        V = V + alpha*PV;
        %r = r-alpha*w
        RU = RU - alpha*WU;
        RV = RV - alpha*WV;
        %计算delta
        delta = sqrt(norm(RU,'fro')^2+norm(RV,'fro')^2);
    end


end











