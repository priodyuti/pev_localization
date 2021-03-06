%%
%% Priodyuti Pradhan, Complex Syatems Lab 
%% priodyutipradhan@gmail.com, 
%% Indian Institute of Technology Indore
function [a_New] = Edge_Rewiring(a, N, Edges)
  T = ((N*(N-1))/2) - Edges;
  k = 1;
  l = 1;
  e1 = zeros(Edges,2);   %Edges present
  e2 = zeros(T,2);       %Edges not present
  for i = 1:(N-1) 
    for j = i+1:N 
      if a(i,j) ==1
        e1(k,1) = i;
        e1(k,2) = j;
        k = k + 1;
      else
        e2(l,1) = i;
        e2(l,2) = j;
        l = l + 1;
      end
    end                               
  end
  %fprintf('k=%d l=%d\n',k,l);
 %e1
 %e2 
  a_New = zeros(N,N);
  while(1)
   a_New = a;
   r1 = randi([1 (k-1)],1,1);
   a_New(e1(r1,1),e1(r1,2)) = 0;  
   a_New(e1(r1,2),e1(r1,1)) = 0;
    
   r0 = randi([1 (l-1)],1,1);
   a_New(e2(r0,1),e2(r0,2)) = 1; 
   a_New(e2(r0,2),e2(r0,1)) = 1;
   
   %-------checking for connectedness-----------
   A = sparse(a_New);
   d = dfs(A,1);
   r = ismember(-1,d); 
   if r ~= 1 
     %fprintf('Network is connected\n'); 
     clear e1 e2
     return 
   end
   %-----------END-----------------------------
 end
end

