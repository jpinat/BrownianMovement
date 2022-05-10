#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define Frec_Max 50
#define Tr_Max 1000000

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

float dist[Tr_Max],dist_final[Tr_Max],frec_dist[Frec_Max];
FILE *f;
extern float Random(void);
extern void ini_ran(int SEMILLA);

main(){
float alfa,x,y,corte,inc,min,max,delta,dist_N;
int NITER,N_tr,i,j,inter;
//alfa=orden intervalo del movimiento
//deltah=distancia recorrida en x o y por paso
//dist_N=distancia tras N pasos
//dist_final=distancia en t final
//frec_dist=histograma para dist_final
//inc=pasito que das al avanzar (delta_h)
//inter=inter del histograma que utilizo
//Norma=normalizo frec_dist
//dist[]=distancia en i al (0,0)

ini_ran(time(0));

alfa=10.0; //Cuanto mayor alfa mejor poisson
corte=100.0; //Utilizo corte como N_intervalos
NITER=500;
N_tr=500000;

for(i=0;i<Frec_Max;i++){ //Inicializo histograma de distancia
frec_dist[i]=0;
}
for(i=0;i<Tr_Max;i++){ //Inicializo de nuevo
dist[i]=dist_final[i]=0;
}
for(j=0;j<N_tr;j++){
x=y=0.0; //Vuelvo al (0,0) para cada nueva trayectoria
for(i=0;i<NITER;i++){
inc=2*alfa*(0.5-Random());
x+=inc;
inc=2*alfa*(0.5-Random());
y+=inc;
dist_N=sqrt(x*x+y*y);
dist[i]+=dist_N;
}
dist_final[j]+=dist_N; //dist_final como suma de todas las distancias
}

//*************************************HISTOGRAMA*********************************
min=99999999;
max=-9999999;
for(j=0;j<Tr_Max;j++){
if(dist_final[j]<min)min=dist_final[j];
if(dist_final[j]>max)max=dist_final[j];
}
delta=(max-min)/Frec_Max;
for(i=0;i<N_tr;i++){
inter=(int)(dist_final[i]-min)/delta; //inter debe ser entero (int)
frec_dist[inter]++;
}
/*
for(i=0;i<Frec_Max;i++){
printf("%f %f\n",min+i*delta,frec_dist[i]/(delta*N_tray);
}
*/
//***********************************************************************************************

//Escribo la raiz de N de la distribucion
f=fopen("rw_d.txt","wt");
for(i=0;i<NITER;i++){
fprintf(f,"%d %f\n",i,dist[i]/N_tr);
}
fclose(f);//Cierro fichero para abrir con el mismo puntero

//Escribo la distribucion poissoniana
f=fopen("frec_dist.txt","wt");
for(i=0;i<Frec_Max;i++){
fprintf(f,"%f %f\n",min+i*delta,frec_dist[i]/(delta*N_tr));
}
fclose(f);
}

float Random(void){
float r;
ig1=ind_ran-24;
ig2=ind_ran-55;
ig3=ind_ran-61;
irr[ind_ran]=irr[ig1]+irr[ig2];
ir1=(irr[ind_ran]^irr[ig3]);
ind_ran++;
r=ir1*NormRANu;
//printf("r=%f\n",r);
return r;
}

void ini_ran(int SEMILLA){
int INI,FACTOR,SUM,i;
srand(SEMILLA);
INI=SEMILLA;
FACTOR=67397;
SUM=7364893;
for(i=0;i<256;i++){
INI=(INI*FACTOR+SUM);
irr[i]=INI;
}
ind_ran=ig1=ig2=ig3=0;
}
