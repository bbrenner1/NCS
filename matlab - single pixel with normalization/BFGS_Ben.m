function [ x_min, value, gradient ] = BFGS_Ben( x0, func )
%Uses BFGS quasi-Newton algorithm to find the minimum negative log
%likelihood added to the noise contribution. Updates x, gradient of f(x), and
%approximate hessian of f(x) until the gradient fulfills the condition that
%(grad'*grad) is less than or equal to the tolerance (an initialized
%variable). Also included are a gradient function and the backtracking line
%search algorithm to find the optimum step-size, alpha. The function being
%optimized accepts a square matrix, however for the purposes of this
%algorithm, an equivalent vector was used. For this reason, a square matrix
%and vector version of x are maintained.

%Initialize
tol=1e-06;
lengthx=length(x0);
B=eye(lengthx^2);
maxit=25;
xvect=zeros(lengthx^2,1);
gradxnew=gradient2(x0,func);
%gradxnew
i=1;
%gradx

%Turning the initial x guess into a vector
p=1;
for n=1:1:lengthx
    for m=1:1:lengthx
     xvect(p,1)=x0(n,m);
     p=p+1;
    end
end

%Initializing x and x-update
xoldv=xvect;

xnewsq=zeros(lengthx);
xoldsq=x0;

%Beginning BFGS loop
while i<maxit
    
    S=-(B^-1)*gradxnew; %Calculating search direction
    
    
    f=feval(func,xoldsq);
    alpha=backtr(func,f,xoldv,gradxnew,S); %Calculating step-size with backtrack algorithm
    
    %alpha=1;
    xnew=xoldv+alpha*S;

    delx=alpha*S;
    p=1;
    %Updating square matrix of x-update with vector version
    for n=1:1:lengthx
        for m=1:1:lengthx
             xnewsq(n,m)=xnew(p,1);
             p=p+1;
        end
    end
    p=1;
        %Updating square matrix of x with vector version
        for n=1:1:lengthx
             for m=1:1:lengthx
              xoldsq(n,m)=xoldv(p,1);
              p=p+1;
             end
        end
    
    gradxold=gradxnew;
    gradxnew=gradient2(xnewsq,func);
    feval(func,xoldsq);

    
    
    %gradxnew
    %Starting conditional update of Hessian
    if abs(gradxnew'*gradxnew)>tol
        delf=gradient2(xoldsq,func);
        Y=gradient2(xnewsq,func)-delf;
        %Y=gradxnew-gradxold;
        %Y;
        A=(Y*Y')./(Y'*delx);
        %A=(Y*Y')*((Y'*delx)^-1);
        %A;
        
        %gradxold;
        %xoldsq;
        %delf;
        C=(delf*delf')./(delf'*S);
        %C=(delf*delf')*((delf'*S)^-1)
       
        B=B+A+C;
       
        i=i+1;
        xoldv=xnew;
        xoldsq=xnewsq;
    else 
        break
    end 
    xoldv=xnew;

%     d=-B*g;
%     alpha=1/((B*s)'*s);
%     xnew=xold+alpha*d;
%     s=xnew-xold;
%     q=feval(func,xnew)-feval(func,xold);
%     B=B+((q*q')/(q'*s))-((B*s*s'*B')/(s'*H*s));
%     %xold=xnew;
%     gradx'*gradx>tol
end
x_min=xnewsq;
value=feval(func,x_min);

gradient=gradxnew;
end
function [ gradx ] = gradient2( xin, func )
hstep = 0.00001;
n = length(xin);
f = feval(func,xin);
p=1;
%gradx=zeros(3);
gradx=zeros(numel(xin),1);
for i1 = 1:n
    for j=1:n
      xs = xin;
      xs(i1,j) = xs(i1,j) + hstep;
      gradx(p,1)= (feval(func,xs) -f)/hstep;
      %gradx(i1,j)= (feval(func,xs) -f)/hstep;
      p=p+1;
    end
end
gradx;
end

function alpha = backtr(objFunc,f,x,gradx,dir)

% 2010 m.bangert@dkfz.de
% backtracking line search using armijo criterion
%
% objFunc      - handle for objective function
% objFuncValue - current objective function value @ x
% x            - x
% dx           - dx
% dir          - search direction
%
% example : mb_backtrackingLineSearch(objFunc,objFuncValue,x,dx,dir)

alphaMax     = 1; % this is the maximum step length
alpha        = alphaMax;
fac          = 1/2; % < 1 reduction factor of alpha
c_1          = 1e-1;
%mu           = 0.2;

lengthx=sqrt(length(x));
xsq=zeros(lengthx);
gradxsq=zeros(lengthx);
p=1;

for n=1:1:sqrt(length(x))
    for m=1:1:sqrt(length(x))
              xsq(n,m)=x(p,1);
              p=p+1;
    end
end
p=1;

for n=1:1:sqrt(length(gradx))
    for m=1:1:sqrt(length(gradx))
        gradxsq(n,m)=gradx(p,1);
        p=p+1;
    end
end
p=1;
dirsq=zeros(sqrt(length(x)));
for n=1:1:sqrt(length(dir))
    for m=1:1:sqrt(length(dir))
        dirsq(n,m)=dir(p,1);
        p=p+1;
    end
end

%if 
%objFunc(xsq+dirsq*alpha)
%objFunc(xsq)
%c_1*alpha*(gradxsq'*dirsq)
    %while abs(objFunc(xsq-c_1*gradxsq)) <=abs(objFunc(xsq)) - abs(c_1*alpha*((sqrt(sum(gradx.^2)))^2))
   
    while objFunc(xsq+dirsq*alpha) > abs(f) + c_1*alpha*(dir'*gradx)
        if (objFunc(xsq+dirsq*alpha)-f) <.000001
            break 
        else 
         alpha = fac*alpha;
        end 
            if alpha < 10*eps
                error('Error in Line search - alpha close to working precision');
            end
    end
% else
%         while objFunc(xsq+dirsq*alpha) <= objFunc(xsq) - c_1*alpha*(gradxsq'*dirsq)
%     
%       alpha = fac*alpha;
%     
%         if alpha < 10*eps
%             error('Error in Line search - alpha close to working precision');
%         end
%     
%     end
% end
end
