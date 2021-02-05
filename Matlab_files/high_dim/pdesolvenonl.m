function [x,grid] = pdesolvenonl(l,y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input: 
    % l = control of the PDE discretization
    % y = unknown parameter
    %%% Output:   
    % x = solution of the PDE
    % grid = grid for the PDE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % determine mesh width
    N=2^l;
    h=1/N;
    x0 = 0;
    x1 = 1;
    grid = x0:h:x1;

    %check size of y
    [m,J] = size(y);
    x = zeros(m,J);

    % rhs.
    rhs=@(x) 1;

    % initialization of A

    for i=1:J
        a=zeros(N-1,1);
        n=zeros(N-2,1);
        co=coeff(y(:,i));
        for k=1:N-2
            tauk1=co((k),1);
            tauk2=co((k+1),1);
            tauk3=co((k+2),1);
            a(k,1)=1/h*(tauk1+2*tauk2+tauk3)/2;
            n(k,1)=-1/h*(tauk2+tauk3)/2;
        end
        tauk1=co((N-1),1);
        tauk2=co((N),1);
        tauk3=co((N+1),1);
        a(N-1,1)=1/h*(tauk1+2*tauk2+tauk3)/2;
        % initialization of rhs f
        f = zeros(1,length(grid)-2);
        for k=1:N-1
            f(k)=rhs(k*h)*h;
        end

        A=diag(a) + diag(n,-1) + diag(n,1);
        % solve the linear system Ax=b
        x(:,i)=A\f';
    end
end

function u=coeff(y)
% coefficient function lognormal distribution
u=exp(y);
u=[0; u; 0];
end




