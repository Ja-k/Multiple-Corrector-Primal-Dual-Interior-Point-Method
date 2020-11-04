
% With Factorization
function [ primal , x , status ] = PDIPM(Struct , MaxIter)

if ~isstruct(Struct)
    error('SQP is not a struct!');
end
if ~isscalar(MaxIter) || ~isreal(MaxIter)
    error('iter is not a real scalar!');
end
if MaxIter <= 0
    error('iter must be a > 0');
end

if ~ isfield( Struct , 'Q' ) || ~ isfield( Struct , 'q' ) || ...
   ~ isfield( Struct , 'E' ) || ~isfield(Struct,'x') || ~isfield(Struct,'b') ||...
 ~isfield(Struct,'lambda')
   error( 'SQP not a well-formed struct' );
end

if ~ isreal( Struct.Q ) || ~ isreal( Struct.q ) || ...
   ~ isreal( Struct.E ) || ~isreal(Struct.x) || ~isreal(Struct.b) ||...
 ~isreal(Struct.lambda)
   error( 'SQP not a well-formed struct' );
end

%output = fopen('output.txt','w');
%--------------------------------------------------------------------------
%    Initializing Parameters
%--------------------------------------------------------------------------

x = Struct.x;
b = Struct.b;
n = numel(x);
k = numel(b);
%x = ones(n,1);
%x = Struct.en;

Q = Struct.Q;
q = Struct.q;

E = Struct.E;
lambda = Struct.lambda;
%lambda = zeros(k,1);



M = zeros( k , k );
e = ones(n,1);
z = Struct.z;
%z = ones(n,1);

% reducing rate of complimentarity value (mu)
epp = 1 / n;

fprintf('Primal-Dual Interior Point Method\n\n');
fprintf('Iter\t\tPrimal\t\tRelative_duality_gap\t\tPrimal_infeas\t\tDual_infeas\t\tComplementarity\n\n');
%fprintf(output,'Iter\t\t\tPrimal\t\t\tDual\tRelative_duality_gap\t\tPrimal_infeas\t\tDual_infeas\n\n');

i = 1;
% Main iteration
while true
    
    xQx = ( x' * Q ) * x;
    
    r_d =  ( Q * x ) + q - ( E' * lambda ) - z;
    
    
    % Primal Objective value
    primal_1 = ( xQx ) + ( q' * x ) ;
    
    primal = 0.5 * ( xQx ) + ( q' * x ) ;
    
    
    % Dual Objective value
    %dual = ( b' * lambda ) - xQx ;  
    
    % Compute primal infeasibility "should be less than 1e-8" for
    % optimality
    primal_infeasibility =  ( norm ( ( E * x ) - b) )  / ( norm ( b ) + 1 );
    
    % Compute dual infeasibilty  "should be less than 1e-6" for optimality
    dual_infeasibility = norm ( r_d ) / ( norm (q) + 1 );
    %dual_infeasibility = (norm ( (E * x) - b , 1 )) / (norm(x,1) + 1);
    %dual_infeasibility = norm( (E * x) - b  ) / ( norm ( ( Q * x ) + q ) + 1 ) ;
        
    % Compute relative dualty gap "should be less than 1e-8" for optimalty
    relative_duality_gap = ( primal_1 - ( b' * lambda  ) ) / ( abs ( primal_1 ) + 1 );
    %relative_duality_gap = abs(((q' * x) + (x' * Q * x )-(b' * lambda)))/(abs((b' * lambda) - (x'* Q * x )) + 1);
    
    %---------------------------------------------------------------------

    % complimentarity value 
    mu = epp * ( ( x' * z ) / n );
    mu_e = mu*e;
    
    

    fprintf( '%4d\t\t%1.8e\t\t%1.4e\t\t%1.4e\t\t%1.4e\t\t%1.4e\n' , i , primal , relative_duality_gap,primal_infeasibility,dual_infeasibility,(x'*z));
   
    % Stop if it reaches MaxIter
    if i > MaxIter
        status = 'Stopped';
        fprintf("x' * z = %d \n", (x' * z));
        break;
    end
    
   % Stop when optimal solution found
    if (x' * z) < 1.000000e-08 
        status = 'optimal';
        break;
    end   
        
    % Diagonal matrice of x
    X = diag(x);
    
    % Diagonal matrice of z
    Z = diag(z);    
    
    
    XZe = sparse(X) * sparse(Z) * e;
    %================================== COEFFICIENT MATRIX ================
    % solve the Linear system "Newton's method" with indefinitefactorization 
    % Forming matrix A of linear equation
    XZ = sparse(X) \ sparse(Z) ;
    T = Q + XZ ;
    
    % Coefficient Matrix
    A = [-T,E';E,M];
    %fprintf("\nCondition of A = %d \n\n", cond(A));
    %fprintf("Condition of A = %d \n", condest(A));
    %fprintf("Determinant of matrix A = %d \n", det(A));
    %--------------- Coefficient matrix Factorization----------------------
    
    [L,D,P] = ldl(A);
%    dA = decomposition(A);
%    tf = isIllConditioned(dA);
%    disp(tf);
  
   %========================== Predictor ==================================
   % Right hand side for predictor 
   %B = [(Q*x)+q-(E'*lambda)-(X\((mu*e))) ; b-(E*x)] ;
   B = [r_d-(X\((mu_e)-XZe));b-(E*x)];
   %B = [(Q*x)+q-(E'*lambda)-z;b-(E*x)];
   %B = [zeros(n,1); zeros(k,1)] ;
   
   
   % Solving The Linear System
   Delta_pred = P * (L' \ (D \ (L \ (P' * B))));
   
%    bb = P' * B;
%    cc = L \ bb;
%    dd = D \ sparse(cc);
%    ee = L' \ dd;
%    Delta_pred = P * ee;
   
   %delta = A \ B;
   %fprintf("norm by back slash = %d \n", norm(Delta_pred - delta) / norm(delta));
   %Delta_pred = delta;
   DeltaX_pred = Delta_pred(1:n,1);
   %DeltaL_pred = Delta_pred(n+1:n+k,1);
   DeltaZ_pred = X \ ( ( mu_e ) - XZe - (sparse(Z) * sparse(diag(DeltaX_pred)) * e) );
  
   %=============================Corrector=================================
   for f = 1:3
       
       dz = X\((mu_e) - XZe - (sparse(diag(DeltaX_pred)) * sparse(diag(DeltaZ_pred)) * e ));
       
        % Right hand side for corrector
        
        %BB = [(Q*x)+q-(E'*lambda)-(X\(mu*e));b-(E*x)];
        BB = [r_d-dz;b-(E*x)];
        %BB = [(Q*x)+q-(E'*lambda)-z;b-(E*x)];
        %BB = [zeros(n,1);zeros(k,1)];
     
        % sovling augmented Linear system
        Delta = P * (L' \ (D \ (L \ ( P' * BB))));
        
%         bb = P' * BB;
%         cc = L \ bb;
%         dd = D \ sparse(cc);
%         ee = L' \ dd;
%         Delta = P * ee;
        
        
        %delta_corr = A\BB;
        %Delta = delta_corr;
        %fprintf("norm by back slash _ corr = %d \n", norm(Delta - delta_corr) / norm(delta_corr));
      
        DeltaX = Delta(1:n,1);
        DeltaL = Delta(n+1:n+k,1);
        DeltaZ = dz - (X\(sparse(Z) * sparse(diag(DeltaX)) * e)) ;     

        DeltaX_pred = DeltaX ;
        DeltaZ_pred = DeltaZ ;
       
    end
    
    %====================================================================== 
       % Calculating the step length
        alpha_primal = ( 1 / ( max ( 1 , max ( - ( DeltaX ./ x ) ) ) ) ) ;
        alpha_dual   = ( 1 / ( max ( 1 , max ( - ( DeltaZ ./ z ) ) ) ) ) ;
    
        
        % Unified step length
        alpha = min ( alpha_primal , alpha_dual );
  
        % Taking steps towards the direction
        x = x + ( 0.99995 * alpha * DeltaX );
        z = z + ( 0.99995 * alpha * DeltaZ ); 
        lambda = lambda + ( 0.99995 * alpha * DeltaL);
       
       
    % check x > 0 and z > 0
    
    if ~all ( x > 0 ) 
        error('negatives in x !');
    end 
     
    if ~all ( z > 0 )
        error('negatives in z !');
    end   
    
    % Increment Iteration
    i = i + 1;
% End of Main iteration    
end  
%fclose(output);
