clear all
load('x1.mat')
load('x2.mat')

n=5;% n>k, n needd large to show correctness
t=2; % maximum is: floor((k)/2);
k=4; % 
no_codewords=n^2; % countmax is to the size of polynomial p(n)


%generate random key x or weight t
rx = [ones(1,t), zeros(1,k-t)]; 
x=rx( randperm(k) )';



y_prime=[];trial=0; 
while isempty(y_prime) && trial<=no_codewords
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Encoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,matrix_F]=encode_T(x1',n,k,x);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[G_prime,x_prime,y_prime]=decoding(x2',n,k,t,y,matrix_F);


if isempty(y_prime) && trial==no_codewords
    disp('no solution')

end 
trial=trial+1;
disp(['trial ', num2str(trial-1), ' of ',num2str((no_codewords))])
end


En=['Encoded Key  : ' int2str(x)'];
Rec=['Recovered Key: ' int2str(x_prime)'];
 disp(En)
 disp(Rec)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G,x,y]=Decode_T(binary_M,k,error_max,query_codeword)
error_weight=error_max;
% d=n-k+1;%define 2^m=n-1
% binary_M=double(binary_M.x);



nstop=nchoosek((k),error_weight);
% binary_M_xor=gf(binary_M);
e_vector_i_2=[zeros(1,(k)-error_weight),ones(1,error_weight)];
jj=1;

while jj<=nstop  % -1 for excluding the zero vector
    e_vector_i_2=e_vector_i_2';
    %     B2test_jj=gf(e_vector_i_2);
    Test_solution=binary_M*e_vector_i_2;
    solution=mod(Test_solution,2);
    
    dist_vector=solution~=query_codeword;
    
    if sum(dist_vector)<=error_weight
        jj=nstop+1;
        
        % disp(['decoding t= ', num2str(error_weight),'/',num2str(error_max)])
        %                 disp('YES, solution found')
        
        x=e_vector_i_2;
        G=binary_M;
        y=solution';
        %uncomment below for global code
        % [G,y,y_query]=get_rows(binary_M_xor,dist_vector,solution,query_codeword);
        
        
    else
        G=[];
        x=[];
        y=[];
        %                     disp(['decoding t=: ', num2str(error_weight),'/',num2str(error_max)])
        
        
    end
    jj=jj+1;
    e_vector_i_2=nextperm(e_vector_i_2);
end


if isempty(G) && isempty(x) && isempty(y) && error_weight==error_max+1
    %     disp('no solution')
end
end
%%
function [ query_codeword,proj_n_mat]=encode_T(x1,n,k,rx)
[binary_M,proj_n_mat]=project_function(x1,n,k); % input ned row vector


query_codeword=binary_M*rx;
query_codeword=mod(query_codeword,2);
end
%%
function [binary_M,proj_n_mat]=project_function(input_Strig,n,k)
proj_n_mat=randn(n*k,size(input_Strig,1));
strin_n_mat=proj_n_mat*input_Strig;
binary_M=reshape(strin_n_mat,n,k);
binary_M=sign(binary_M);
binary_M(binary_M==-1)=0;
end
%%
function [binary_M]=project_function2(input_Strig,n,k,mat1)
strin_n_mat=mat1*input_Strig;
binary_M=reshape(strin_n_mat,n,k);
binary_M=sign(binary_M);
binary_M(binary_M==-1)=0;
end
%%
function [G,x,y]=decoding(x2,n,k,t,query_codeword,proj_n_mat)
[binary_M2]=project_function2(x2,n,k,proj_n_mat); % input ned row vector
% binary_M2=gf(binary_M2);
[G,x,y]=Decode_T(binary_M2,k,t,query_codeword);
end
%%
