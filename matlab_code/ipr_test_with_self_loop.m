%%
clear
clc

bins1 = 24;
bins2 = 40;
path ='/home/priodyuti/Desktop/self-loop/ER-100/Networks/ER_1.txt';

fd = fopen(path,'rt');

formatSpec = '%d %d';
sizeA = [2 Inf];

Y = fscanf(fd, formatSpec, sizeA);
Y = Y';
fclose(fd);

N = max(max(Y(:,1)),max(Y(:,2)));
A = zeros(N, N);
t = size(Y, 1);

for i=1:t
    A(Y(i,1),Y(i,2)) = 1;
    A(Y(i,2),Y(i,1)) = 1;
end    
clear Y;
%A=A+eye(N,N);
deg = sum(A,2);
[max_deg max_index] = max(deg);
total_deg = sum(A(:));

fprintf('Number of Nodes: %d\n', N);
fprintf('Number of edges: %d\n',total_deg/2 ); 
fprintf('Max_Deg: %d\n', max_deg); 
fprintf('Max_Index: %d\n', max_index); 
fprintf('Average_deg: %d\n', total_deg/N); 
%fprintf('Determinant of A: %d\n', det(A)); 
%fprintf('Condition number of A: %d\n', cond(A));

[evec, eval] = eig(A);
cv = diag(eval);

fprintf('Max_Eigenval: %f\n', cv(N));
fprintf('Second largest eigenval: %f\n', cv(N-1));

Ev1 = evec(:, N);

IPR = sum(Ev1.^4);
fprintf('IPR = %f\n', IPR);

%L = diag(sum(A)) - A;
L = A + diag(sum(A));
%D = diag(sum(A));
%d = diag(D);
%for i=1:N
%  D(i,i) = 1/d(i);
%end
%T = A*D;

%[evec1, eval1] = eig(A);
[evec2, eval2] = eig(L);
%[evec3, eval3] = eig(T);
cv2 = diag(eval2);

d = sqrt(diag(L));
R = L./(d*d'); %correlation matrix

[evec3, eval3] = eig(R);
cv3 = diag(eval3)


%C=corrcov(L)
%[ sigma, corr ] = cov2corr(L);


%clear all


%cv

%norm(A*evec(:,1)-eval(1)*evec(:,1))

%figure(1)
%plot(1:N,cv,'r.');
% sort(diag(eval),'descend');

% eigen_vector1=sort(abs(Ev),'descend');
% for i=1:N
%  fprintf('%0.32f\n',abs(Ev(i)));
% end

%eigen_vec=sort(abs(Ev),'descend');
%figure(2)
% plot(1:N,eigen_vec,'r.');
%hist(eigen_vec,50)
%clear A


