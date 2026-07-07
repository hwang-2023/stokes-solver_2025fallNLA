function [U2, V2]= GS(U, V, F, G, N)%实现GS迭代

    h = 1/N;
   
    %更新速度分量U
    [U1,V1] = left_multiply_A(U,V,N); 


    RF=F-U1;

    EU = zeros(N+1,N);

    EU(1,1) = RF(1,1);
    for i=2:N
        EU(i,1) = 1/3*EU(i-1,1) + 1/3*h^2*RF(i,1);
    end
    EU(N+1,1) = RF(N+1,1);

    for j = 2:N-1
        EU(1,j) = RF(1,j);
        for i = 2:N
            EU(i,j) = 1/4*EU(i,j-1) + 1/4*EU(i-1,j) + 1/4*h^2*RF(i,j);
        end
        EU(N+1,j) = RF(N+1,j);
    end
    
    EU(1,N) = RF(1,N);
    for i = 2:N
        EU(i,N) = 1/3*EU(i-1,N) + 1/3*EU(i,N-1) + 1/3*h^2*RF(i,N);
    end
    EU(N+1,1) = RF(N+1,1);

    U2 = U + EU;
    


    RG=G-V1;
    
    
    EV = zeros(N,N+1);


    EV(1,1) = RG(1,1);
    for j=2:N
        EV(1,j) = 1/3*EV(1,j-1) + 1/3*h^2*RG(1,j);
    end
    EV(1,N+1) = RG(1,N+1);

    
    
    for i=2:N-1
        EV(i,1) = RG(i,1);
        for j=2:N
            EV(i,j) = 1/4*EV(i,j-1) + 1/4*EV(i-1,j) + 1/4*h^2*RG(i,j);
        end
        EV(i,N+1) = RG(i,N+1);
    end

    EV(N,1) = RG(N,1);
    for j=2:N
        EV(N,j) = 1/3*EV(N,j-1) + 1/3*h^2*RG(N,j);
    end
    EV(N,N+1) = RG(N,N+1);

    
    
    V2 = V + EV;




    
end




