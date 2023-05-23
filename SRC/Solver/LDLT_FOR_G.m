%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     Perform LDLT for Gene_Alpha_Integ_Linear                    *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Gene_Alpha_Integ_Linear.m                                 *
%*                                                                 *
%*                                                                 *
%* - Programmed by:                                                *
%*     Yang Mo                                                     *
%*     Institute of Systems Engineering                            *
%*     China Academy of Engineering Physic, 2023.05.12             *
%*                                                                 *
%* *****************************************************************
function[a]=LDLT_FOR_G(m,q)
% 函数为对ma=q做LDLT分解并求解向量a
% m为待LDLT分解的矩阵
% q为列向量
%%%%%%对矩阵M进行LDLT分解%%%%%%%%
n=size(m,1);
a=zeros(n,1);    %%%解向量      
L = eye(n);      %%%L矩阵
u=zeros(n);
D = zeros(n,1);  %%%因为D矩阵为对角阵，因此用n行列向量进行存储，节省存储空间
D(1)=m(1,1);
for i=2:n
    for j=1:i-1
        temp=0;
        for k=1:j-1
            temp=temp+u(i,k)*L(j,k);
        end
        u(i,j)=m(i,j)-temp;
        L(i,j)=u(i,j)/D(j);
    end
    temp=0;
    for k=1:i-1
        temp=temp+u(i,k)*L(i,k);
    end
    D(i)=m(i,i)-temp;
end
%D=diag(D);
%%%%%%对矩阵M进行LDLT分解%%%%%%%%
%%%%%%%%%%%%用LDLT分解求解线性方程组m*a=q%%%%%%%%%%%%%%%%%%%%%
temp1=zeros(n,1);     
temp3=zeros(n,1);
for j=1:n
    temp2=0;
    for k=1:j-1
        temp2=temp2+L(j,k)*temp1(k);
    end
    temp1(j)=q(j)-temp2;
    temp3(j)=(q(j)-temp2)/D(j);
end
    for j=1:n
        temp4=0;
        for k=n-j+2:n
            temp4=temp4+L(k,n-j+1)*a(k);
        end
        a(n-j+1)=temp3(n-j+1)-temp4;  
    end
%%%%%%%%%%%%用LDLT分解求解线性方程组m*a=q%%%%%%%%%%%%%%%%%%%%%
end