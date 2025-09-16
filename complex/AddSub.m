%%

% create two complex numbers
a = complex(4,5);
b = 3+2i;

% let MATLAB do the hard work
z1 = a+b;

% the "manual" way
z2 = complex( real(a)+real(b) , imag(a)+imag(b) );

%% subtraction is the same as addition...

% let MATLAB do the hard work
z3 = a-b;

% the "manual" way
z4 = complex( real(a)-real(b) , imag(a)-imag(b) );

%% done.
