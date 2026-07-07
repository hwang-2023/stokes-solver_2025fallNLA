function [U, V, iter, time] = Uzawa(F, G, n, alpha, MAXITER)

    %进行Uzawa迭代
    %MAXITER:最大迭代次数

    tic;

    N = 2^n;
    r0 = sqrt(norm(F,'fro')^2 + norm(G,'fro')^2);

    U = zeros(N+1,N);
    V = zeros(N,N+1);
    P = zeros(N,N);

    rh = r0;
    
    
    iter = 0;

    while rh/r0 > 1e-8 && iter < MAXITER
        iter = iter + 1;
        
        %AU = F - BP
        [F1,G1] = left_multiply_B(P,N);
        [U,V,iter1] = CG(F-F1,G-G1,N,1e-8);

        %P = P + alpha * (BT*U)
        P = P + alpha*left_multiply_BT(U,V,N);
  
        [U1,V1,P1] = left_multiply_matrix(U,V,P,N);
        rh = sqrt(norm(U1-F,'fro')^2 + norm(V1-G,'fro')^2 + norm(P1,'fro')^2);

        disp(['Uzawa iteration',num2str(iter),': CG iter=',num2str(iter1),', rh=',num2str(rh)]);

    end


    time = toc;
end
