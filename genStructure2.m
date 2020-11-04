% getStructure function for creating the instances with two parameters
% from the user
function Struct = genStructure2( n , k)

if ~isscalar(n) || ~isreal(n)
    error(' n is not a real scalar !');
end
if ~isscalar(k) || ~isreal(k)
    error(' k is not a real scalar !');
end
if n <= 0 || k <= 0
    error('n and k must be greater than 0 ');
end
if k > n 
    error('k must be < n');
end    

% Random sparse matrix for sparse instaces
G = sprandn( n - k , n , 0.01);

% Random dense matrix for dense instaces
%G = randn( n - k  , n);

% Semdefinite nxn matrix for objective function from randomly created matrices
Q = G' * G ;
Struct.Q = Q ;


% Random nx1 vector as q vector in objective function
q = randn( n , 1 );
Struct.q = q;


% kx1 matrix of ones for linear constraint of objective function
b = ones( k , 1);
Struct.b = b;

% nx1 matrix of ones to be used in computations
e = ones(n,1);


%------------------------------Creating E matrix---------------------------
%--------------------------------------------------------------------------
% Initial kxn matrix of zeros to be modified by following 
U = zeros( k , n );
% counter variable for keep track of modifications
counter = 1;
% assigning k value to a variable to prevent changing the actual value of k
k_value = k;

% assigning n value to a variable to prevent changint the actual value of n
n_value = n ;

% Loop k times (#rows)
for c = 1:k_value
    % in each row the #1's to be computed
    row = ceil(n_value / k_value);
    % modifying U matrix by changing the element value (0) to 1
    for u = 1: row
        U(c,counter) = 1;
        counter = counter + 1;
    end
    % Go to next set of columns
    n_value = n_value - row;
    % Go to next row
    k_value = k_value - 1;
end
Struct.E = U; 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%---------------------------- Computing the initial point------------------
%--------------------------------------------------------------------------
x = U' * ((U * U') \  b) ;
%Struct.x = x;
lambda = ((U * U') \ U) * (Q * x + q);
Struct.lambda = lambda;

z = (Q * x) + q - (U'*lambda) ;

dx = max(-(1.5)* min(x),0);
dz = max(-(1.5) * min(z),0);

x = x + dx * e;
z = z + dz * e;

ddx = 0.5 * ((x' * z)/(e' * z));
ddz = 0.5 * ((x' * z)/(e' * x));

x = x + ddx;
z = z + ddz;

Struct.x = x;
Struct.z = z;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

