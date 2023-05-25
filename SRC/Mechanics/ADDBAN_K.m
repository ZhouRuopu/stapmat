function ADDBAN_K(K, N)

global sdata;
global cdata;

STIFF = sdata.STIFF; IJ = sdata.IJ;
I = IJ(N, 1); J = IJ(N, 2); 

STIFF(6*I-5:6*I, 6*I-5:6*I) = STIFF(6*I-5:6*I, 6*I-5:6*I) + K(1:6, 1:6);
STIFF(6*I-5:6*I, 6*J-5:6*J) = STIFF(6*I-5:6*I, 6*J-5:6*J) + K(1:6, 7:12);
STIFF(6*J-5:6*J, 6*I-5:6*I) = STIFF(6*J-5:6*J, 6*I-5:6*I) + K(7:12, 1:6);
STIFF(6*J-5:6*J, 6*J-5:6*J) = STIFF(6*J-5:6*J, 6*J-5:6*J) + K(7:12, 7:12);

sdata.STIFF = STIFF;