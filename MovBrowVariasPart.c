#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Frec_Max 50 //Numero de intervalos de histograma frec_dist
#define Tr_Max 100000 //Numero maximo de datos a pasar a dicho histograma
#define NormRANu (2.3283063671E-10F)

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

extern double Random(void);
extern void ini_ran(int SEMILLA);

double Random(void){
double r;
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
}

void Histogram(double *data,double *Hist, int N_datos, int N_intervalos,double *del,double *m,double *M){
int i,Indice;
double Norma,max,min,delta;
//1ºInicializo el vector Hist(i<N_intervalos, sino Hist[] tiene inicializadas mas componentes a 0 de las que tiene y luego Indice se repite todo el rato
for(i=0;i<N_intervalos;i++)
    Hist[i]=0;
//2ºCalculo el maximo y minimo de los numeros generados
max=*M;
min=*m;
delta=*del;
for(i=0;i<N_datos;i++){
    if(data[i]>max) max=data[i];
    if(data[i]<min) min=data[i];
}
//Calculo delta, el histograma y la norma
delta=(max-min)/N_intervalos;
for(i=0;i<N_datos;i++){
        Indice=(int)(data[i]-min)/delta;
        Hist[Indice]++;
}
Norma=1.0/(N_datos*delta);
*M=max;
*m=min;
*del=delta;
//Normalizo el vector Hist
for(i=0;i<N_intervalos;i++){
    Hist[i]=Hist[i]*Norma;
}
}



main(){
double alfa,corte,deltah,frec_dist[Frec_Max],dist[Tr_Max],dist_final[Tr_Max],dist_N,min,max,delta,Norma;
//alfa=orden intervalo del movimiento
//deltah=distancia recorrida en x o y por paso
//dist_N=distancia tras N pasos
//dist_final=distancia en t final
//frec_dist=histograma para dist_final
int i,N_pas,N_tray,j;
ini_ran(time (0));
FILE*f1,*f2,*f3;
f1=fopen("rw_d.txt","w");
f2=fopen("frec_dist.txt","w");
f3=fopen("movimiento.txt","w");
alfa=1.0;
corte=100.0;
N_pas=500; //pasos que doy por trayectoria
N_tray=20; //numero de trayectorias que hago
double x[N_tray],y[N_tray];
min=999999999;
max=-9999999999;
for(i=0;i<Frec_Max;i++){
    frec_dist[i]=0;
}
for(i=0;i<N_tray;i++){
    dist[i]=dist_final[i]=0;
    x[i]=y[i]=0;
}
//*****************************************************
for(j=0;j<N_pas;j++){

    for(i=0;i<N_tray;i++){
        deltah=2*alfa*(0.5-Random()); //Genero un deltah y lo sumo en x
        x[i]=x[i]+deltah;
        deltah=2*alfa*(0.5-Random()); //Hago lo mismo para y, generando un deltah distinto al de x
        y[i]=y[i]+deltah;
        dist_N=sqrt(x[i]*x[i]+y[i]*y[i]);
        dist[i]+=dist_N;
        fprintf(f3,"%lf %lf\t",x[i],y[i]);
    }
        fprintf(f3,"\n");
        dist_final[j]+=dist_N;
}
Histogram(dist_final,frec_dist,N_tray,Frec_Max,&delta,&min,&max);
printf("delta=%lf\n",delta);
Norma=1.0/(delta*N_tray);
for(i=0;i<Frec_Max;i++){
    //printf("frec_dist[%d]=%lf\n",i,frec_dist[i]);
}
for(i=0;i<N_pas;i++){
fprintf(f1,"%d %lf\n",i,dist[i]);
}
for(i=0;i<Frec_Max;i++){
fprintf(f2,"%lf %lf\n",min+i*delta,frec_dist[i]);
}
    fclose(f1);
    fclose(f2);
    fclose(f3);
}

