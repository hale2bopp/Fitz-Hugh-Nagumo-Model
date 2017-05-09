function dydt = oscnetwork_opt_vec_out_fast(t,y)
    % declare global variables
    global N g I E c 
    
    % set constants
    vth = 0.3; vsl = 0.001;
    A = 3; B= 3; e = 0.01; c = -0.1;
    
    % create an empty vector to hold the system of differential equations
    dydt = zeros(N*(N+1),1);
    
    % create a matrix holding the synaptic variables
    synvar =reshape(y(2*N+1:N*(N+1)),N-1,N);
   
    % matrix with diagonal as reversal potentials
    excvar = eye(N);
    excvar(1:N+1:N^2) = E;
    
    % insert a zero diagonal into the synaptic variable matrix
    synvarB = [zeros(1,N); tril(synvar)] + [triu(synvar); zeros(1,N)];
    synvarB(1:N+1:N^2)=0;
    
    % matrix with diagonal as fast variables
    Y =eye(N);
    Y(1:N+1:N^2) =y(1:N);
    
    % find the sum for the fast variable

    sum1=diag(g*(synvarB*Y)) - diag(g*(synvarB'*excvar)');
    
    % set fast variables
%    dydt(1:N) = (y(1:N).*(y(1:N)+c).*(1-y(1:N)) - y(N+1:2*N) + ((0.5+0.5*sign(5-t))*I)-sum1)/e;
    dydt(1:N) = (y(1:N).*(y(1:N)+c).*(1-y(1:N)) - y(N+1:2*N) + (I)-sum1)/e;
    % set recovery variables
    dydt(N+1:2*N) = y(1:N) - 0.5.*y(N+1:2*N);

    % set synaptic variables
    for ii =1:N    
        dydt((ii-1)*(N-1)+2*N+1:(ii-1)*(N-1)+2*N+1+N-2) =  (A*  0.5*( (1+ atan((y(1:N~=ii)-vth)/vsl) ) ).*(1-synvar(:,ii)))-B*synvar(:,ii);            
    end
end
