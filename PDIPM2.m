function [ primal , x , status ] = PDIPM2 (Struct , MaxIter)

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
   ~ isreal( Struct.E ) || ~ isreal(Struct.x) || ~isreal(Struct.b) ||...
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
x = ones(n,1);
%x = Struct.en;

Q = Struct.Q;
q = Struct.q;

E = Struct.E;
%lambda = Struct.lambda;
lambda = zeros(k,1);



M = zeros( k , k );
e = ones(n,1);
%z = Struct.z;
z = ones(n,1);
epp = 0.5;

fprintf('Primal-Dual Interior Point Method\n\n');
fprintf('Iter\t\t\tPrimal\t\t\tDual\tRelative_duality_gap\t\tPrimal_infeas\t\tDual_infeas\n\n');
%fprintf(output,'Iter\t\t\tPrimal\t\t\tDual\tRelative_duality_gap\t\tPrimal_infeas\t\tDual_infeas\n\n');

i = 1;
% Main iteration
while true
    
    % Diagonal matrice of x
    %X = diag(x);
    
    % Diagonal matrice of z
    %Z = diag(z);
    
    % Primal Objective value
    primal = (((x' * Q) * x)/2)+ (q' * x);
    
    % Dual Objective value
    dual = (b' * lambda)-((x' * Q) * x);  
    
    % Compute primal infeasibility "should be less than 1e-6" for
    % optimality
    primal_infeasibility = ( norm ( ( E * x ) - b ) ) / ( norm ( b ) + 1 );
    
    % Compute dual infeasibilty  "should be less than 1e-6" for optimality
    dual_infeasibility = norm( q + ( Q * x ) - ( E' * lambda )-z ) / (norm(q )+1);
    %dual_infeasibility = norm( (E * x) - b ) / (norm( (Q * x )+ q )+1);
    %dual_infeasibility = (primal - dual)/ abs(dual);
    
    % Compute relative dualty gap "should be less than 1e-8" for optimalty
    relative_duality_gap = ((q' * x) + ((x' *Q) * x )-(b' * lambda))/(abs((q' * x)+ ((x'* Q) * x )) + 1 );
    %relative_duality_gap = ((q' * x) + (x' * Q * x )-(b' * lambda))/(abs((b' * lambda) - (x'* Q * x )) + 1);
    
    %---------------------------------------------------------------------
%     rp = b - (E * x);
%     rd = (Q * x)+q - (E' * lambda)-z;
   % dgap = (x' * Q * x)+(q' * x)-(b' * lambda);
   % disp(dgap);
    %---------------------------------------------------------------------
   % rpp = norm(rp,1)/(norm(x,1)+1);
   % disp(rpp);
   % rdd = norm(rd,1)/(norm(z,1)+1);
   % disp(rdd);
   % cc = ((x'*z)/n)/(norm(x,1)+ norm(z,1)+1);
   % disp(cc);
    %tol = 1.000000e-8;
%     if (rpp < tol) && (rdd < tol) && (cc<tol)
%         status = 'optimal by tol';
%         break;
%     end
    if (x' * z) < 1.000000e-8
        status = 'optimal by comp.';
        break;
    end     
    mu = epp * ( ( x' * z ) / n );
   % disp(mu);
%     if mu < 1.000000e-08
%         status = 'mu';
%         break;
%     end
    fprintf( '%4d\t\t%1.8e\t\t%1.8e\t\t%1.4e\t\t%1.4e\t\t%1.4e\n' , i , primal , dual , relative_duality_gap,primal_infeasibility,dual_infeasibility);
   % fprintf(output,'%4d\t\t%1.8e\t\t%1.8e\t\t%1.4e\t\t%1.4e\t\t%1.4e\n' , i , primal , dual , relative_duality_gap,primal_infeasibility,dual_infeasibility);
    if i > MaxIter
        status = 'Stopped';
        break;
    end
    % check wether infeasibilities and duality gap satity the conditions
   % if ((primal_infeasibility < 1.000000e-6) && (dual_infeasibility < 1.000000e-06)) && (relative_duality_gap < 1.000000e-08)
    %if (relative_duality_gap < 1.000000e-08)
    %if mu < 1.000000e-08
      %  fprintf(output,'\nThe optimal point is: \n');
      %  fprintf(output,'%1.4f \n', x);
      %  fprintf(output,'\n');
    %    status = 'Optimal';
     %   break;
    %else
        
    % Diagonal matrice of x
    X = diag(x);
    
    % Diagonal matrice of z
    Z = diag(z);    
    %================================== COEFFICIENT MATRIX ===============
    % solve the Linear system "Newton's method" with indefinitefactorization 
    % Forming matrix A of linear equation
    ll = -(Q + ( X \ Z));
    %A( 1:n , 1:n ) = -(Q + ( X \ Z));
    A( 1:n , 1:n ) = ll;
    A(1:n , (n+1):(n+k) ) = E';
    A( (n+1):(n+k) , 1:n )= E ;
    A( (n+1):(n+k) , (n+1):(n+k))= M ;
    
    %--------------- Coefficient matrix Factorization----------------------
    [L,D,P] = ldl(A);
   
    
   %==========================Affine   mu = 0 =============================
   
   B_affine(1:n,1)=(Q * x ) + q - (E' * lambda) ;
   B_affine((n+1):(n+k),1) = b - (E * x) ;
   bb = P' * B_affine;
   cc = L\bb;
   dd = D\cc;
   ee = L'\dd;
   Delta_affine = P * ee;
   DeltaX_affine=Delta_affine(1:n,1);
   %DeltaL_affine=Delta_affine((n+1):(n+k) , 1 );
   DeltaZ_affine= -z - ((X \ Z) * DeltaX_affine) ;
   
   %========================Affine step length=============================
   % affine primal step length with upperbound = 1
   alphax_affine = (1/(max ( 1 , max((-DeltaX_affine ./ x))))) ;
   
   % affine dual step length with upperbound = 1
   alphaz_affine=  (1/(max ( 1 , max((-DeltaZ_affine ./ z)))));
   
   % mu_affine computation
   mu_affine = ((x + (alphax_affine * DeltaX_affine))' * (z + (alphaz_affine * DeltaZ_affine))) / n;
   
   % computing cetering parameter
   epp = (mu_affine / mu)^3;
   %==========================Solving Step Direction system================
    B(1:n,1) = (Q * x ) + q  - (E' * lambda)-z ;
    B((n+1):(n+k),1) = b - (E * x) ;
    %tic
    bb = P'*B;
    cc = L\bb;
    dd = D\cc;
    ee = L'\dd;
    Delta = P*ee;
    %toc
    DeltaX = Delta( 1:n , 1 );
    DeltaL = Delta( ( n+1 ):( n+k ) , 1 );
    DeltaZ = ( mu * ( X \ e) )- z - ( X \ (Z * DeltaX )) - (X\(DeltaX_affine .* DeltaZ_affine)) ;
    
    %================================== step length========================
    % compute the step length
    alpha_primal = 0.99995 * ( 1 / (max ( 1 , max(-(DeltaX ./ x) ) ) ) ) ;
    alpha_dual = 0.99995 * (1/(max ( 1 , max ( - ( DeltaZ ./ z) ) ) ) );
    alpha = min(alpha_primal , alpha_dual );
    
    % update the vector x, z, lambda
     x = x + ( alpha * DeltaX );
     z = z + ( alpha * DeltaZ );
     lambda = lambda + ( alpha * DeltaL );
    %======================================================================  
    % check x > 0 and z > 0
    if ~ all(x > 0) 
         error('negatives in x !');
     end    
    if ~ all(z > 0 )
        error('negatives in z !');
    end   
    
    %end
    
    % Increment Iteration
    i = i + 1;
% End of Main iteration    
end  
%fclose(output);
