for (boot = 1; boot < 101; boot += 1) {
/*Import data*/

DataSet myData         = ReadDataFile("/bootstap_data/data_" + boot + ".fasta");
DataSetFilter myFilter = CreateFilter(myData,3,"","","TGG,TGT,TTA,TTC,TTG,TTT");
HarvestFrequencies (F, myFilter, 3, 3, 1);


/*Define global parameters*/

global Paux1 = 0.5; Paux1 :< 1;
global Paux2 = 0.5; Paux2 :< 1;
global Paux3 = 0.5; Paux3 :< 1;

global P0 := Paux1 * Paux2;
global P1 := (1-Paux1) * Paux3;
global P2 := (1-Paux1) * (1-Paux3);
global P3 :=  Paux1 * (1-Paux2);

global m10=0.00000000000000000001; m10 :< 10; 
global m12=0.00000000000000000001; m12 :< 10;  
global m13=0.00000000000000000001; m13 :< 10;  

global m20=0.00000000000000000001; m20 :< 10;  
global m21=0.00000000000000000001; m21 :< 10;  
global m23=0.00000000000000000001; m23 :< 10;  

global m30=0.00000000000000000001; m30 :< 10;  
global m31=0.00000000000000000001; m31 :< 10;  
global m32=0.00000000000000000001; m32 :< 10;

global s_del=-0.1; s_del :> - 1;

global s0 := (1-2*s_del);
global s1 := (1-s_del);
global s2 := 1;
global s3 := (1-(s_del/2));
  

m01 :=0;  
m02 :=0;
m03 :=0;
m13 :< m12;
m20 :< m21;
m30 :< m31;


Ppolvar:=1.0;

/*Define ancestral frequencies*/
Ppol := ((F[4]+F[5]+F[6]+F[7]+F[8]+F[9]+F[10]+F[11]+F[12]+F[13]+F[14]+F[15]+F[16]+F[17]+F[18]+F[19]+F[20]+F[21]+F[22]+F[23]+F[24]+F[25]+F[26]+F[27]+F[28]+F[29]+F[30]+F[31]+F[32]+F[33]+F[34]+F[35]+F[36]+F[37]+F[38]+F[39]+F[40]+F[41]+F[42]+F[43]+F[44]+F[45]+F[46]+F[47]+F[48]+F[49]+F[50]+F[51]+F[52]+F[53]+F[54]+F[55]+F[56]+F[57]) * Ppolvar);

norm := ( P0 * (m01 + m02 + m03) ) + ( P3 * (m31 + m32 + m30) ) + ( P1 * (m10 + m12 + m13) ) + ( P2 * (m20 + m21 + m23) );

Freqs={{P0*(1-Ppol)}
{P1*(1-Ppol)}
{P2*(1-Ppol)}
{P3*(1-Ppol)}
{Ppol*((P0*m01*0.111111111111111) + (P1*m10*1))/norm}
{Ppol*((P0*m01*0.125) + (P1*m10*0.5))/norm}
{Ppol*((P0*m01*0.142857142857143) + (P1*m10*0.333333333333333))/norm}
{Ppol*((P0*m01*0.166666666666667) + (P1*m10*0.25))/norm}
{Ppol*((P0*m01*0.2) + (P1*m10*0.2))/norm}
{Ppol*((P0*m01*0.25) + (P1*m10*0.166666666666667))/norm}
{Ppol*((P0*m01*0.333333333333333) + (P1*m10*0.142857142857143))/norm}
{Ppol*((P0*m01*0.5) + (P1*m10*0.125))/norm}
{Ppol*((P0*m01*1) + (P1*m10*0.111111111111111))/norm}
{Ppol*((P0*m02*0.111111111111111) + (P2*m20*1))/norm}
{Ppol*((P0*m02*0.125) + (P2*m20*0.5))/norm}
{Ppol*((P0*m02*0.142857142857143) + (P2*m20*0.333333333333333))/norm}
{Ppol*((P0*m02*0.166666666666667) + (P2*m20*0.25))/norm}
{Ppol*((P0*m02*0.2) + (P2*m20*0.2))/norm}
{Ppol*((P0*m02*0.25) + (P2*m20*0.166666666666667))/norm}
{Ppol*((P0*m02*0.333333333333333) + (P2*m20*0.142857142857143))/norm}
{Ppol*((P0*m02*0.5) + (P2*m20*0.125))/norm}
{Ppol*((P0*m02*1) + (P2*m20*0.111111111111111))/norm}
{Ppol*((P0*m03*0.111111111111111) + (P3*m30*1))/norm}
{Ppol*((P0*m03*0.125) + (P3*m30*0.5))/norm}
{Ppol*((P0*m03*0.142857142857143) + (P3*m30*0.333333333333333))/norm}
{Ppol*((P0*m03*0.166666666666667) + (P3*m30*0.25))/norm}
{Ppol*((P0*m03*0.2) + (P3*m30*0.2))/norm}
{Ppol*((P0*m03*0.25) + (P3*m30*0.166666666666667))/norm}
{Ppol*((P0*m03*0.333333333333333) + (P3*m30*0.142857142857143))/norm}
{Ppol*((P0*m03*0.5) + (P3*m30*0.125))/norm}
{Ppol*((P0*m03*1) + (P3*m30*0.111111111111111))/norm}
{Ppol*((P1*m12*0.111111111111111) + (P2*m21*1))/norm}
{Ppol*((P1*m12*0.125) + (P2*m21*0.5))/norm}
{Ppol*((P1*m12*0.142857142857143) + (P2*m21*0.333333333333333))/norm}
{Ppol*((P1*m12*0.166666666666667) + (P2*m21*0.25))/norm}
{Ppol*((P1*m12*0.2) + (P2*m21*0.2))/norm}
{Ppol*((P1*m12*0.25) + (P2*m21*0.166666666666667))/norm}
{Ppol*((P1*m12*0.333333333333333) + (P2*m21*0.142857142857143))/norm}
{Ppol*((P1*m12*0.5) + (P2*m21*0.125))/norm}
{Ppol*((P1*m12*1) + (P2*m21*0.111111111111111))/norm}
{Ppol*((P1*m13*0.111111111111111) + (P3*m31*1))/norm}
{Ppol*((P1*m13*0.125) + (P3*m31*0.5))/norm}
{Ppol*((P1*m13*0.142857142857143) + (P3*m31*0.333333333333333))/norm}
{Ppol*((P1*m13*0.166666666666667) + (P3*m31*0.25))/norm}
{Ppol*((P1*m13*0.2) + (P3*m31*0.2))/norm}
{Ppol*((P1*m13*0.25) + (P3*m31*0.166666666666667))/norm}
{Ppol*((P1*m13*0.333333333333333) + (P3*m31*0.142857142857143))/norm}
{Ppol*((P1*m13*0.5) + (P3*m31*0.125))/norm}
{Ppol*((P1*m13*1) + (P3*m31*0.111111111111111))/norm}
{Ppol*((P2*m23*0.111111111111111) + (P3*m32*1))/norm}
{Ppol*((P2*m23*0.125) + (P3*m32*0.5))/norm}
{Ppol*((P2*m23*0.142857142857143) + (P3*m32*0.333333333333333))/norm}
{Ppol*((P2*m23*0.166666666666667) + (P3*m32*0.25))/norm}
{Ppol*((P2*m23*0.2) + (P3*m32*0.2))/norm}
{Ppol*((P2*m23*0.25) + (P3*m32*0.166666666666667))/norm}
{Ppol*((P2*m23*0.333333333333333) + (P3*m32*0.142857142857143))/norm}
{Ppol*((P2*m23*0.5) + (P3*m32*0.125))/norm}
{Ppol*((P2*m23*1) + (P3*m32*0.111111111111111))/norm}};


/*Define matrices: each rate depends on global mutational parameters (mac etc.), on the global selection coefficient sc (and the synonym sg), and on the branch-specific parameter t representing branch lengths (see HyPhy manual for this last one). From these, population genetic parameters are estimated as explained in the main text.*/ 

matrix1 = {{*,0,0,0,0,0,0,0,0,0,0,0,t*m01,0,0,0,0,0,0,0,0,t*m02,0,0,0,0,0,0,0,0,t*m03,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,*,0,0,t*m10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*m12,0,0,0,0,0,0,0,0,t*m13,0,0,0,0,0,0,0,0,0}
{0,0,*,0,0,0,0,0,0,0,0,0,0,t*m20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*m21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*m23}
{0,0,0,*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*m30,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*m31,0,0,0,0,0,0,0,0,t*m32,0,0,0,0,0,0,0,0}
{0,t*10*1*9/((((1+(s0-s1))*1) + 9)*10),0,0,*,t*10*1*9*(1+(s0-s1))/((((1+(s0-s1))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,t*10*2*8/((((1+(s0-s1))*2) + 8)*10),*,t*10*2*8*(1+(s0-s1))/((((1+(s0-s1))*2) + 8)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,t*10*3*7/((((1+(s0-s1))*3) + 7)*10),*,t*10*3*7*(1+(s0-s1))/((((1+(s0-s1))*3) + 7)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,t*10*4*6/((((1+(s0-s1))*4) + 6)*10),*,t*10*4*6*(1+(s0-s1))/((((1+(s0-s1))*4) + 6)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,t*10*5*5/((((1+(s0-s1))*5) + 5)*10),*,t*10*5*5*(1+(s0-s1))/((((1+(s0-s1))*5) + 5)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s0-s1))*6) + 4)*10),*,t*10*6*4*(1+(s0-s1))/((((1+(s0-s1))*6) + 4)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s0-s1))*7) + 3)*10),*,t*10*7*3*(1+(s0-s1))/((((1+(s0-s1))*7) + 3)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s0-s1))*8) + 2)*10),*,t*10*8*2*(1+(s0-s1))/((((1+(s0-s1))*8) + 2)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{t*10*9*1*(1+(s0-s1))/((((1+(s0-s1))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s0-s1))*9) + 1)*10),*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,t*10*1*9/((((1+(s0-s2))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,*,t*10*1*9*(1+(s0-s2))/((((1+(s0-s2))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*2*8/((((1+(s0-s2))*2) + 8)*10),*,t*10*2*8*(1+(s0-s2))/((((1+(s0-s2))*2) + 8)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*3*7/((((1+(s0-s2))*3) + 7)*10),*,t*10*3*7*(1+(s0-s2))/((((1+(s0-s2))*3) + 7)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*4*6/((((1+(s0-s2))*4) + 6)*10),*,t*10*4*6*(1+(s0-s2))/((((1+(s0-s2))*4) + 6)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*5*5/((((1+(s0-s2))*5) + 5)*10),*,t*10*5*5*(1+(s0-s2))/((((1+(s0-s2))*5) + 5)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s0-s2))*6) + 4)*10),*,t*10*6*4*(1+(s0-s2))/((((1+(s0-s2))*6) + 4)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s0-s2))*7) + 3)*10),*,t*10*7*3*(1+(s0-s2))/((((1+(s0-s2))*7) + 3)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s0-s2))*8) + 2)*10),*,t*10*8*2*(1+(s0-s2))/((((1+(s0-s2))*8) + 2)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{t*10*9*1*(1+(s0-s2))/((((1+(s0-s2))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s0-s2))*9) + 1)*10),*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,t*10*1*9/((((1+(s0-s3))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,*,t*10*1*9*(1+(s0-s3))/((((1+(s0-s3))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*2*8/((((1+(s0-s3))*2) + 8)*10),*,t*10*2*8*(1+(s0-s3))/((((1+(s0-s3))*2) + 8)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*3*7/((((1+(s0-s3))*3) + 7)*10),*,t*10*3*7*(1+(s0-s3))/((((1+(s0-s3))*3) + 7)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*4*6/((((1+(s0-s3))*4) + 6)*10),*,t*10*4*6*(1+(s0-s3))/((((1+(s0-s3))*4) + 6)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*5*5/((((1+(s0-s3))*5) + 5)*10),*,t*10*5*5*(1+(s0-s3))/((((1+(s0-s3))*5) + 5)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s0-s3))*6) + 4)*10),*,t*10*6*4*(1+(s0-s3))/((((1+(s0-s3))*6) + 4)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s0-s3))*7) + 3)*10),*,t*10*7*3*(1+(s0-s3))/((((1+(s0-s3))*7) + 3)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s0-s3))*8) + 2)*10),*,t*10*8*2*(1+(s0-s3))/((((1+(s0-s3))*8) + 2)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{t*10*9*1*(1+(s0-s3))/((((1+(s0-s3))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s0-s3))*9) + 1)*10),*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,t*10*1*9/((((1+(s1-s2))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,*,t*10*1*9*(1+(s1-s2))/((((1+(s1-s2))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*2*8/((((1+(s1-s2))*2) + 8)*10),*,t*10*2*8*(1+(s1-s2))/((((1+(s1-s2))*2) + 8)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*3*7/((((1+(s1-s2))*3) + 7)*10),*,t*10*3*7*(1+(s1-s2))/((((1+(s1-s2))*3) + 7)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*4*6/((((1+(s1-s2))*4) + 6)*10),*,t*10*4*6*(1+(s1-s2))/((((1+(s1-s2))*4) + 6)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*5*5/((((1+(s1-s2))*5) + 5)*10),*,t*10*5*5*(1+(s1-s2))/((((1+(s1-s2))*5) + 5)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s1-s2))*6) + 4)*10),*,t*10*6*4*(1+(s1-s2))/((((1+(s1-s2))*6) + 4)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s1-s2))*7) + 3)*10),*,t*10*7*3*(1+(s1-s2))/((((1+(s1-s2))*7) + 3)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s1-s2))*8) + 2)*10),*,t*10*8*2*(1+(s1-s2))/((((1+(s1-s2))*8) + 2)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,t*10*9*1*(1+(s1-s2))/((((1+(s1-s2))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s1-s2))*9) + 1)*10),*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,t*10*1*9/((((1+(s1-s3))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,*,t*10*1*9*(1+(s1-s3))/((((1+(s1-s3))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*2*8/((((1+(s1-s3))*2) + 8)*10),*,t*10*2*8*(1+(s1-s3))/((((1+(s1-s3))*2) + 8)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*3*7/((((1+(s1-s3))*3) + 7)*10),*,t*10*3*7*(1+(s1-s3))/((((1+(s1-s3))*3) + 7)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*4*6/((((1+(s1-s3))*4) + 6)*10),*,t*10*4*6*(1+(s1-s3))/((((1+(s1-s3))*4) + 6)*10),0,0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*5*5/((((1+(s1-s3))*5) + 5)*10),*,t*10*5*5*(1+(s1-s3))/((((1+(s1-s3))*5) + 5)*10),0,0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s1-s3))*6) + 4)*10),*,t*10*6*4*(1+(s1-s3))/((((1+(s1-s3))*6) + 4)*10),0,0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s1-s3))*7) + 3)*10),*,t*10*7*3*(1+(s1-s3))/((((1+(s1-s3))*7) + 3)*10),0,0,0,0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s1-s3))*8) + 2)*10),*,t*10*8*2*(1+(s1-s3))/((((1+(s1-s3))*8) + 2)*10),0,0,0,0,0,0,0,0,0}
{0,t*10*9*1*(1+(s1-s3))/((((1+(s1-s3))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s1-s3))*9) + 1)*10),*,0,0,0,0,0,0,0,0,0}
{0,0,0,t*10*1*9/((((1+(s2-s3))*1) + 9)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,*,t*10*1*9*(1+(s2-s3))/((((1+(s2-s3))*1) + 9)*10),0,0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*2*8/((((1+(s2-s3))*2) + 8)*10),*,t*10*2*8*(1+(s2-s3))/((((1+(s2-s3))*2) + 8)*10),0,0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*3*7/((((1+(s2-s3))*3) + 7)*10),*,t*10*3*7*(1+(s2-s3))/((((1+(s2-s3))*3) + 7)*10),0,0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*4*6/((((1+(s2-s3))*4) + 6)*10),*,t*10*4*6*(1+(s2-s3))/((((1+(s2-s3))*4) + 6)*10),0,0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*5*5/((((1+(s2-s3))*5) + 5)*10),*,t*10*5*5*(1+(s2-s3))/((((1+(s2-s3))*5) + 5)*10),0,0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*6*4/((((1+(s2-s3))*6) + 4)*10),*,t*10*6*4*(1+(s2-s3))/((((1+(s2-s3))*6) + 4)*10),0,0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*7*3/((((1+(s2-s3))*7) + 3)*10),*,t*10*7*3*(1+(s2-s3))/((((1+(s2-s3))*7) + 3)*10),0}
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*8*2/((((1+(s2-s3))*8) + 2)*10),*,t*10*8*2*(1+(s2-s3))/((((1+(s2-s3))*8) + 2)*10)}
{0,0,t*10*9*1*(1+(s2-s3))/((((1+(s2-s3))*9) + 1)*10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,t*10*9*1/((((1+(s2-s3))*9) + 1)*10),*}};



dim = Rows (matrix1);
sparse_matrix = {dim,dim};
Model M1 = (matrix1, Freqs, 0);

for (r = 0; r < dim; r += 1) {
    for (c = 0; c < dim; c+=1) {
        if (r != c) {
            GetString (info, M1, r, c); 
            if (Abs (info) > 0 && info != "0") {
                ExecuteCommands ("sparse_matrix[r][c]:=" + info);
            }
        }
    }
}


/*Defining tree and constraints on tree*/

Model M_sparse = (sparse_matrix, Freqs,0);


ACCEPT_ROOTED_TREES=1;
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;
Tree myTree=((pop_1:1.0,pop_3:1.0)a1:1.0,(pop_2:1.0,pop_4:1.0)a2:1.0);

myTree.a1.t :> 0.1; myTree.a2.t :> 0.1; myTree.pop_1.t :> 0.05; myTree.pop_2.t :> 0.05;



/*Likelihood optimization*/

LikelihoodFunction3 Fun = (myFilter, myTree, Freqs); 
Optimize(results,Fun); 


/*printing the output*/

fprintf(stdout,"After optimization: ", Fun);

out="/HyPhy_output.txt";

fprintf(out, "\n|", boot ,"|\n");
fprintf(out, Fun);


}

