function  f= gridSearch(func,b,s,data,mat,show)

%%
%b = points on the grid
%s = 1 x size(b) matrix giving the number of points tried for each range
%Performs a linear grid search



gridSize=max(s);
bgridValues=zeros(gridSize,size(b,1));

for i=1:size(b,1)
   temp=zeros(gridSize-s(i),1);
   bgridValues(:,i)=[linspace(b(i,1),b(i,2),s(i))';temp];
end




sizeb=size(bgridValues,2);
loop=prod(s(1:sizeb-mat));
intresults=zeros(loop,sizeb+1);
rowsS=length(s);

if mat>0
   gmat=createGrid(bgridValues(:,(sizeb-mat+1):sizeb),s((sizeb-mat+1):sizeb));
end

if sizeb-mat>0
   
if mat>0
for i=1:loop
   	gtemp=[bgridValues(1+floor((i-1)/(prod(s(1:rowsS-mat))/s(1))),1)*ones(prod(s(sizeb-mat+1:sizeb)),1)];
      if sizeb-mat>1
         for j=2:sizeb-mat
            t=mod(1+floor((i-1)/(prod(s(j:rowsS-mat))/s(j))),s(j));
            if t==0
               t=s(j);
            end
            gtemp=[gtemp, bgridValues(t,j)*ones(prod(s(sizeb-mat+1:sizeb)),1)];   
         end
      end
   gtemp=[gtemp, gmat];
 	out=feval(func, gtemp',data);
   temp1=sortrows([out, gtemp],1);
	intresults(i,:)= temp1(1,:);  
   clear out temp1 gtemp;
   
    if(show)
	    disp(strcat('loop iteration  : ',num2str(i), ', of  :', num2str(loop)));
    end
end  %Loop

else   %%mat<=0
   for i=1:loop
   	gtemp=[bgridValues(1+floor((i-1)/(prod(s(1:rowsS))/s(1))))];
      if sizeb>1
         for j=2:sizeb-mat
            t=mod(1+floor((i-1)/(prod(s(j:rowsS))/s(j))),s(j));
            if t==0
               t=s(j);
            end
            gtemp=[gtemp, bgridValues(t,j)];   
         end
      end
 	out=feval(func, gtemp',data);
   temp1=sortrows([out, gtemp],1);
	intresults(i,:)= temp1(1,:);  
   clear out temp1 gtemp;
	if(show)
	    disp(strcat('loop iteration  : ',num2str(i), ', of  :', num2str(loop)));
    end

	end %Loop
end  %%mat>0



else		%%%sizeb-mat==0
   out=feval(func,gmat',data);
   intresults=[out,gmat];
   clear out gtemp;
end		


f.results=sortrows(intresults,1);
f.val=f.results(1,1);
f.b=f.results(1,2:size(f.results,2))';





















