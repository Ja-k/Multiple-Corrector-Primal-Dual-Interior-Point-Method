function [DeltaX , DeltaZ] = predictor(L,P,D,B,X,Z,XZe,mu_e,e,n)

    Delta_pred = P * (L' \ (D \ (L \ (P' * B))));
    DeltaX = Delta_pred(1:n,1);
    DeltaZ = X \ ( ( mu_e ) - XZe - (Z * sparse(diag(DeltaX)) * e) );
  