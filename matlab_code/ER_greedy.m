%%
%% Priodyuti Pradhan, Complex Syatems Lab 
%% priodyutipradhan@gmail.com, 
%% Indian Institute of Technology Indore
function  ER_greedy(N, k, Iter)
  clc
  fprintf('Number of nodes: %d\n', N);
  fprintf('Average degree: %d\n', k);
  fprintf('Number of Iterations: %d\n', Iter);
  p = k/N;
  fprintf('Connection probability: %f\n', p);
  %A = zeros(N,N);

  path1 = 'ipr.txt'; 
  path2 = 'ER_opt.txt'; 

  fd1 = fopen(path1,'wt');
  fd2 = fopen(path2,'wt'); 


  [A Edges] = ER_Random(N, p);
  fprintf('Number of edges: %d\n\n', Edges);
  %A
  [evec1,eval1] = eig(A);
  evals = diag(eval1);
  max_eig_val = evals(N);
  sec_max_eig_val = evals(N-1);
  c1 = sum(evec1(:,N).^4);
  %[max_deg max_index] = max(sum(A,2));

  clear evec1 eval1
  Flag = 0;
  fprintf(fd1,'%d %f %f %f %f %f\n',Flag, c1, c1, max_eig_val, sec_max_eig_val, min(evals));
 
  for iter = 2:Iter 
    %fprintf('before\n');
    Flag = 0;
    [A_New] = Edge_Rewiring(A, N, Edges);
  %  A_New

    cvx_begin sdp
      variable d(N)
      minimize(sum(d))
      subject to
         A_New + diag(d) >= 0
         d >= 0
    cvx_end

    A_New = A_New + diag(d); %covariance matrix

    x = sqrt(diag(A_New));
    R = A_New./(x*x'); %correlation matrix  


    [evec1 eval1] = eig(R);
    evals = diag(eval1);
    max_eig_val = evals(N);
    sec_max_eig_val = evals(N-1);
    c2 = sum(evec1(:,N).^4);

    clear evec1 eval1 R
    
    if c2 > c1
      A = A_New; 
      %for i=1:N
      %   A(i,i) = 0;
      %end
      c1 = c2;
      Flag = 1; 
    end
  
    fprintf(fd1,'%d %f %f %f %f %f\n', Flag, c1, c2, max_eig_val, sec_max_eig_val, min(evals));      
    clear A_New
  end
  %A
  %hist(diag(A))
  for i=1:N
    for j=i:N
      if A(i,j) ~= 0
        fprintf(fd2,'%d %d %f\n',i,j,A(i,j));
      end 
    end
  end
  fclose(fd1);
  fclose(fd2);
  clear all
  fprintf('\nCompletes...\n');
 
end
