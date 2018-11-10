%%
%% Priodyuti Pradhan, Complex Syatems Lab 
%% priodyutipradhan@gmail.com, 
%% Indian Institute of Technology Indore
function [g edges] = ER_Random(N, p)
  g = zeros(N,N); 
  while(1)
   edges = 0;
   for i=1:N-1
     for j=i+1:N;
        x = rand;  
        %cout<<x<<endl; 
        if x<p
          g(i,j) = 1;
          g(j,i) = 1;
          edges = edges + 1;
        else
          g(i,j) = 0;
          g(j,i) = 0;
        end
     end
     g(i,i) = 0;
   end
   %------Checking for connectedness-------
   A = sparse(g);
   d = dfs(A,1);
   r = ismember(-1,d);
   if r ~= 1  %nu == N
     %fprintf('Network is connected\n');
     return
   end 
   %--------- END--------------------------
  end
end







