function [DeltaX , DeltaZ] = predictor(L,P,D,B_p,X,Z,XZe,mu_e,e,n)

    Delta_pred = P * (L' \ (D \ (L \ (P' * B_p))));
    DeltaX = Delta_pred(1:n,1);
    DeltaZ = X \ ( ( mu_e ) - XZe - (Z * sparse(diag(DeltaX)) * e) );
  