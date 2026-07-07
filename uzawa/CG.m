function [U, V, iter] = CG(F, G, N, epsilon)
    
    %(实用)共轭梯度法求解AU = F - BP
    %F,G:RHS的矩阵形式
    %epsilon:停止条件


    b = sqrt(norm(F,'fro')^2+norm(G,'fro')^2);
    rho = b^2;
    iter = 0;

    U = zeros(N+1,N);
    V = zeros(N,N+1);

    RU = F;
    RV = G;

    %取一个最大共轭梯度法迭代次数
    MAXITER = N*sqrt(N);

    while sqrt(rho) > epsilon*b && iter < MAXITER 
        iter = iter + 1;
        if iter == 1
            %k=1, p=r
            PU = RU;
            PV = RV;
            rho = norm(RU,'fro')^2 + norm(RV,'fro')^2;
        else
            %k>1时, 计算beta, 再求下降方向p=r+beta*p
            rho_before = rho;
            rho = norm(RU,'fro')^2 + norm(RV,'fro')^2;
            beta = rho/rho_before;
            PU = RU + beta*PU;
            PV = RV + beta*PV;
        end
        %w=Ap
        [WU,WV] = left_multiply_A(PU,PV,N);
        alpha = rho/(dot_multiply(WU,PU)+dot_multiply(WV,PV));
        %x=x+alpha*p
        U = U + alpha*PU;
        V = V + alpha*PV;
        %r=r-alpha*w
        RU = RU - alpha*WU;
        RV = RV - alpha*WV;
        
    end


end











