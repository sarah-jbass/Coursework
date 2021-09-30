function V=numVarianceCov(func,b,data,W,e)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given a function "func" that calculates the objective function
%and moment conditions, this function calculates the se of the
%parameters of the function.

index=[1:length(b)]';

f=feval(func,b,data);

for i=1:length(b)
   u=index==(i*ones(length(b),1));
   fH=feval(func,b+b(i)*e*u,data);
   M(:,i)=((fH.m-f.m)./(b(i)*e));
end

omega=f.omega;

V.VarCov=(1/f.N)*inv(M'*W*M)*M'*W*omega*W*M*inv(M'*W*M);
V.se=diag(V.VarCov).^(.5);





