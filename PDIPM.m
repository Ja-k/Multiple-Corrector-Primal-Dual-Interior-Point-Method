
% Primal Dual Interior Points Method
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

%--------------------------------------------------------------------------
%    Initializing Parameters
%--------------------------------------------------------------------------

x = Struct.x;
b = Struct.b;
n = numel(x);
k = numel(b);


Q = Struct.Q;
q = Struct.q;

E = Struct.E;
lambda = Struct.lambda;

M = zeros( k , k );
e = ones(n,1);
z = Struct.z;

epp = 1/n;

fprintf('Primal-Dual Interior Point Method\n\n');
fprintf('Iter\t\tPrimal\t\tRelative_duality_gap\t\tPrimal_infeas\t\tDual_infeas\t\tComplementarity\n\n');
%corr_time = 0 ;

i = 1;
% Main iteration
while true
    
    xQx = ( x' * Q ) * x;
    
    r_d =  ( Q * x ) + q - ( E' * lambda ) - z;
    
    
    % Primal Objective value
    primal_1 =( xQx ) + ( q' * x ) ;
    
    primal = 0.5 * ( xQx ) + ( q' * x ) ;
    
    
    % Dual Objective value
    %dual = ( b' * lambda ) - xQx ;  
    
    % Compute primal infeasibility "should be less than 1e-8" for
    % optimality
    primal_infeasibility =  ( norm ( ( E * x ) - b) )  / ( norm ( b ) + 1 );
    
    % Compute dual infeasibilty  "should be less than 1e-6" for optimality
    dual_infeasibility = norm ( r_d ) / ( norm (q) + 1 );
        
    % Compute relative dualty gap "should be less than 1e-8" for optimalty
    relative_duality_gap = ( primal_1 - ( b' * lambda  ) ) / ( abs ( primal_1 ) + 1 );
    
    %---------------------------------------------------------------------

    % complimentarity value 
    mu = epp * ( ( x' * z ) / n );
    
    mu_e = mu * e;
    
    

    fprintf( '%4d\t\t%1.8e\t\t%1.4e\t\t%1.4e\t\t%1.4e\t\t%1.4e\n' , i , primal , relative_duality_gap,primal_infeasibility,dual_infeasibility,(x'*z));
   
    % Stop if it reaches MaxIter
    if i > MaxIter
        status = 'Stopped';
        fprintf("x' * z = %d \n", (x' * z));
        break;
    end
    
    %roh = max(1,norm(Q),norm(E),norm(q),norm(b));
    %fprintf("r_d value = %d \n",abs(r_d));
   % Stop when optimal solution found
    if ((x' * z) < 1.000000e-8) %&& (abs(r_d) < 1.000000e-8)
        status = 'optimal';
        break;
    end   
        
    % Diagonal matrice of x
    X = sparse(diag(x));
    
    % Diagonal matrice of z
    Z = sparse(diag(z));    
    
    
    XZe = X * Z * e;
    %================================== COEFFICIENT MATRIX ================
    % solve the Linear system "Newton's method" with indefinitefactorization 
    % Forming matrix A of linear equation
    XZ = X \ Z ;
    T = Q + XZ ;
    
    % Coefficient Matrix
    A = [-T,E';E,M];
    %fprintf("Condition of A = %d",cond(A));
    %--------------- Coefficient matrix Factorization----------------------
    
    [L,D,P] = ldl(A);
  
   %========================== Predictor ==================================
   % Right hand side for predictor 
   
   B = [r_d-(X\((mu_e)-XZe));b-(E*x)];   
   
   [DeltaX_pred , DeltaZ_pred] = predictor(L,P,D,B,X,Z,XZe,mu_e,e,n);
   %tic;
   % Solving The Linear System
%    Delta_pred = P * (L' \ (D \ (L \ (P' * B))));
%    DeltaX_pred = Delta_pred(1:n,1);
%    DeltaZ_pred = X \ ( ( mu_e ) - XZe - (Z * sparse(diag(DeltaX_pred)) * e) );
%   
   %toc;
   %=============================Corrector=================================
   %tic;
   
   for f = 1:2
       
%        dz = X\((mu_e) - XZe - (sparse(diag(DeltaX_pred)) * sparse(diag(DeltaZ_pred)) * e ));
%        
%         % Right hand side for corrector
%         BB = [r_d-dz ; b - ( E * x ) ];
%         
%         % sovling augmented Linear system
%         Delta = P * (L' \ (D \ (L \ ( P' * BB))));
%       
%         DeltaX = Delta(1:n,1);
%         DeltaL = Delta(n+1:n+k,1);
%         DeltaZ = dz - (X\(Z * sparse(diag(DeltaX)) * e)) ;     

        
        [DeltaX, DeltaL , DeltaZ] = corrector(L,P,D,X,Z,XZe,mu_e,e,n,E,x,DeltaX_pred,DeltaZ_pred,r_d,b,k);
        
        DeltaX_pred = DeltaX ;
        DeltaZ_pred = DeltaZ ;
       
    end
    %toc;
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
%fprintf("Corr_time = %s \n",sum(corr_time));
