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
    error('k must be <= n');
end    

%G = sprandn( n - k , n , 0.01);
G = randn( n - 10  , n);
%disp(G);
%spy(G);
Q = G' * G;
%spy(Q);
Struct.Q = Q; 
% Generate matrix q
q = randn( n , 1 );
Struct.q = q;
% generate b vector
b = ones( k , 1);
Struct.b = b;

% e vector 
e = ones(n,1);

% generate Lambda vector with only positive values
%--------------------Constructing E matrix according to simplices----------
U = zeros( k , n );
%xx = zeros( n , 1 );

%------------------------------Creating E matrix-------------------
counter = 1;
k_value = k;
n_value = n ;
for c = 1:k_value
    row = ceil(n_value / k_value);
    for u = 1: row
        U(c,counter) = 1;
        counter = counter + 1;
    end
    n_value = n_value - row;
    k_value = k_value - 1;
end
Struct.E = U; 
%disp(U);
%--------------------------------------------------------------------------
% computring initial x
%x = U' * ((U*U')^(-1)) *  b ;
x = U' * ((U*U') \  b) ;

%disp(x);
%Struct.x = x;
%lambda = ((U * U')^(-1)) * U * (Q * x + q);
lambda = ((U * U') \ U) * (Q * x + q);
%disp(lambda);
Struct.lambda = lambda;

z = (Q * x) + q - (U'*lambda) ;
%disp(z);
dx = max(-(1.5)* min(x),0);
dz = max(-(1.5) * min(z),0);

x = x + dx * e;
z = z + dz * e;

ddx = 0.5 * ((x' * z)/(e' * z));
ddz = 0.5 * ((x' * z)/(e' * x));

x = x + ddx;
z = z + ddz;

Struct.x = x;
%disp(Struct.x);
Struct.z = z;
%disp(Struct.z);
