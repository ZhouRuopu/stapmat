function ADDBAN_M(M, N)

global sdata;
global cdata;

MASS = sdata.MASS; IJ = sdata.IJ;
I = IJ(N, 1); J = IJ(N, 2); 

% 行和法化为集中质量阵
M_l = zeros(12,12);
for m = 1:12
    M_l(m,m) = sum(M(1,:));
end

for i = 1:6
MASS(6*I-6+i, 6*I-6+i) = MASS(6*I-6+i, 6*I-6+i) + M_l(i,i);
MASS(6*J-6+i, 6*J-6+i) = MASS(6*J-6+i, 6*J-6+i) + M_l(6+i,6+i);
end

sdata.MASS = MASS;
