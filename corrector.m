function [DeltaX, DeltaL , DeltaZ] = corrector(L,P,D,X,Z,XZe,mu_e,e,n,E,x,DeltaX_pred,DeltaZ_pred,r_d,b,k)

        dz = X\((mu_e) - XZe - (sparse(diag(DeltaX_pred)) * sparse(diag(DeltaZ_pred)) * e ));
       
        % Right hand side for corrector
        BB = [r_d-dz;b-(E*x)];
        
        % sovling augmented Linear system
        Delta = P * (L' \ (D \ (L \ ( P' * BB))));
      
        DeltaX = Delta(1:n,1);
        DeltaL = Delta(n+1:n+k,1);
        DeltaZ = dz - (X\(Z * sparse(diag(DeltaX)) * e)) ;     