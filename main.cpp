#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>

#define H0 -0.5 // Energía inicial, en teoría debería conservarse
#define phi1 0  // ángulo phi*
#define e 0.6 //excentricidad


//funciones


float f_Kepler(float a, float b){  //Ecuacion diferencial de la aceleración 2.2
float f, g;
g=pow((a*a+b*b),1.5);
f= -a/g;
return f;
}

float velocidad(float tiempo, float a, float b){
float v;
v= tiempo*f_Kepler(a,b);
return v;
}

float Energia(float x,float y,float vx,float vy){
    float E;
    E= 0.5*(vx*vx + vy*vy) - 1/(sqrt(x*x + y*y));

return E;
}
void EulerExplicito(int pasos, float h);
void ImpMidpoint(int pasos, float h);
void Verlet(int pasos, float h);
void SympcEuler(int pasos, float h);
void printEexplicito(float x[], float y[],float ener[], float err[], int pasos);
void printEsymp(float x[], float y[],float ener[],float err[], int pasos);
void printMidpoint(float x[], float y[],float ener[],float err[], int pasos);
void printVerlet(float x[], float y[],float ener[],float err[], int pasos);
float Calcerror(float x, float y); ///El valor de d hay que ponerlo a mano
void CalculoAnalitico(int pasos);





using namespace std;

int main()
{
float L0, d;



d= 1.0-e*e;      //definición de d y L0
L0 = sqrt(d);

printf("Parametros del sistema: H0 = %f  e = %f L0 = %f d = %f\n",H0,e, L0, d);
 puts("Pulsa Intro para continuar ... ");
    getchar();

// crear array de los pasos temporales para distintos métodos

float h[4];
int steps[4];
h[0]= 0.0005;
steps[0]=400000 ;//euler explicito
h[1]=0.05;
steps[1]=4000; //symplectic Euler, implicit midpoint, Verlet


///SOLUCION ANALITICA (sale)
//CalculoAnalitico(steps[1]);
puts("Calculo analitico acabado ... ");

///EULER EXPLÍCITO (sale)
//EulerExplicito(steps[0], h[0]);
 puts("Metodo Euler explicito acabado ... ");
 //getchar();

///VERLET (sale)
//Verlet(steps[1],h[1]);
puts("Metodo de Verlet acabado ... ");
getchar();

///IMPLICIT MIDPOINT
ImpMidpoint(steps[1], h[1]);
puts("Metodo implicit midpoint acabado ... ");
getchar();

///SYMPLECTIC EULER (sale)
//SympcEuler(steps[1],h[1]);
puts("Metodo symplectic Euler acabado ... ");
getchar();


    return 0;
}


void EulerExplicito(int pasos, float h){

float x[pasos];
float y[pasos];
float vx[pasos];
float vy[pasos];
float energ[pasos];
float err[pasos];

for (int i=0; i<pasos; i++){  //inicializar vectores
   x[i]=0;
   y[i]=0;
   vx[i]=0;
   vy[i]=0;
   energ[i]=0;
   err[i]=0;
}

x[0]= 1-e;
y[0]=0;
vx[0]=0;
vy[0]=sqrt((1.0 +e)/(1.0 -e));

energ[0]=Energia(x[0],y[0],vx[0],vy[0]);

x[1]= x[0] + h*vx[0];
y[1]= y[0] + h*vy[0];

printf("ha entrado");


for(int i=1; i<pasos-1; i++){
vx[i]= vx[i-1] + h*f_Kepler(x[i-1],y[i-1]); //ojo con el orden al usar f_Kepler
vy[i]=vy[i-1] + h*f_Kepler(y[i-1],x[i-1]);

x[i+1]= x[i] + h*vx[i];
y[i+1]= y[i] + h*vy[i];

energ[i]=Energia(x[i],y[i],vx[i],vy[i]);
err[i]=Calcerror(x[i],y[i]);

printf("Paso %d. X=%f  Y= %f Energia=%f\n", i, x[i], y[i], energ[i]);

}

printEexplicito(x,y,energ,err,pasos);
}

void Verlet(int pasos, float h){

float x[pasos];
float y[pasos];
float vx[pasos];
float vxm[pasos];
float vy[pasos];
float vym[pasos];
float energ[pasos];
float err[pasos];

for (int i=1; i<pasos; i++){  //inicializar vectores
   vx[i]=0;
   vxm[i]=0;
   vy[i]=0;
   vym[i]=0;
   x[i]=0;
   y[i]=0;
   energ[i]=0;
   err[i]=0;

}

x[0]= 1-e;
y[0]=0;
vx[0]=0;
vy[0]=sqrt((1.0 +e)/(1.0 -e));
energ[0]=-0.5;

printf("Valor de h %f \n", h);
//dadas las condiciones iniciales, hay que calcular vm en el punto inicial, para asi poder calcular
//x_i+1 (o y_i+1), pero a la vez, para cañcular vm necesitamos usar h/2 y los puntos x_0, y_0.

for(int i =0; i<(pasos-1);i++){
vxm[i]  = vx[i]+h*0.5*f_Kepler(x[i],y[i]);
vym[i]  = vy[i]+h*0.5*f_Kepler(y[i],x[i]);


x[i+1]= x[i] + h*vxm[i];
y[i+1]=y[i] + h*vym[i];

vx[i+1]= vxm[i] + 0.5*h*f_Kepler(x[i+1],y[i+1]);
vy[i+1]= vym[i] + 0.5*h*f_Kepler(y[i+1],x[i+1]);

//printf("Paso %d. vmx_i = %f --> x_i+1 = %f \n vmy_i = %f --> y_i+1 = %f \n xm_i = %f --> vx_i+1= %f \n ym_i = %f --> vy_i+1 = %f \n", i, vxm[i],x[i+1], vym[i], y[i+1], xm[i], vx[i+1], ym[i], vy[i+1]);
//getchar();
energ[i+1]=Energia(x[i+1],y[i+1],vx[i+1],vy[i+1]);
err[i+1]=Calcerror(x[i+1],y[i+1]);
printf("Paso %d \n", i);

}

printVerlet(x,y,energ,err,pasos);
}



void ImpMidpoint(int pasos, float h){
float x[pasos];
float k1_x[pasos];
float kx[pasos];
float ky[pasos];
float k2_x[pasos];
float y[pasos];
float k1_y[pasos];
float k2_y[pasos];
float vx[pasos];
float k1_vx[pasos];
float k2_vx[pasos];
float vy[pasos];
float k1_vy[pasos];
float k2_vy[pasos];
float energ[pasos];
float err[pasos];
float t=0;
float aux;

for (int i=1; i<pasos; i++){  //inicializar vectores
   vx[i]=0;
   kx[i]=0;
   ky[i]=0;
   k1_vx[i]=0;
   k2_vx[i]=0;
   vy[i]=0;
   k1_vy[i]=0;
   k2_vy[i]=0;
   x[i]=0;
   k1_x[i]=0;
   k2_x[i]=0;
   y[i]=0;
   k1_y[i]=0;
   k2_y[i]=0;
   energ[i]=0;
   err[i]=0;

}

x[0]= 1-e;
y[0]=0.01;
vx[0]=0;
vy[0]=sqrt((1.0 +e)/(1.0 -e));
energ[0]=-0.5;
aux=0;
for(int i=0; i<pasos-1;i++){

k1_vx[i]=f_Kepler(x[i],y[i]);
k1_vy[i]=f_Kepler(y[i],x[i]);

k1_x[i]= x[i] + h*(vx[i] + 0.5*h*k1_vx[i]);
k1_y[i]= y[i] + h*(vy[i] + 0.5*h*k1_vy[i]);

k2_vx[i]=f_Kepler(k1_x[i],k1_y[i]);
k2_vy[i]=f_Kepler(k1_y[i],k1_x[i]);

x[i+1]=x[i] + h*vx[i] + 0.5*h*h*k1_vx[i];
y[i+1]=y[i] + h*vy[i] + 0.5*h*h*k2_vy[i];

vx[i+1]=vx[i]+ 0.5*h*(k1_vx[i]+k2_vx[i]);
vy[i+1]=vy[i]+ 0.5*h*(k1_vy[i]+k2_vy[i]);
/*


k2_vx[i-1]=h*f_Kepler(x[i-1]+0.5*h*vx[i-1], y[i-1]+0.5*h*vy[i-1]);
k2_vy[i-1]=h*f_Kepler(y[i-1]+0.5*h*vy[i-1], x[i-1]+0.5*h*vx[i-1]);

vx[i]=vx[i-1]+k2_vx[i-1];
vy[i]=vy[i-1]+k2_vy[i-1];


k2_x[i-1]=h*(velocidad(h*0.5,x[i-1]+0.5*h*vx[i-1],y[i-1]+0.5*h*vy[i-1]) + vx[i-1]);
k2_y[i-1]=h*(velocidad(h*0.5,y[i-1]+0.5*h*vy[i-1],x[i-1]+0.5*h*vx[i-1]) + vy[i-1]);

x[i]=x[i-1]+k2_x[i-1];
y[i]=y[i-1]+k2_y[i-1];

t=i*h;*/
//printf("tiempo %f",t);
/*
kx[i-1]= x[i-1]+h*0.5*vx[i-1];
ky[i-1]= y[i-1]+h*0.5*vy[i-1];

vy[i]=vy[i-1]+ h*(f_Kepler(ky[i-1],kx[i-1]));
vx[i]=vx[i-1]+ h*(f_Kepler(kx[i-1],ky[i-1]));

k1_vx[i-1]=vx[i-1]+0.5*h*f_Kepler(x[i-1],y[i-1]);
k1_vy[i-1]=vy[i-1]+0.5*h*f_Kepler(y[i-1],x[i-1]);


x[i]=x[i-1] + h*0.5*(vx[i-1]+vx[i]);
y[i]=y[i-1] + h*0.5*(vy[i-1]+vy[i]);*/




energ[i]=Energia(x[i],y[i],vx[i],vy[i]);
err[i]=Calcerror(x[i],y[i]);
//printf("Paso %d. X=%f  Y= %f Energia=%f  INI X0 = %f Y0 =%f\n", i, x[i], y[i], energ[i], x[0],y[0]);
//getchar();
}

printMidpoint(x,y,energ,err,pasos);
}

//imprimir archivos de texto con las coordenadas x,y en los distintos pasos de tiempo, y también la energía
void SympcEuler(int pasos, float h){

float x[pasos];
float y[pasos];
float vx[pasos];
float vy[pasos];
float energ[pasos];
float err[pasos];

for (int i=0; i<pasos; i++){  //inicializar vectores
   x[i]=0;
   y[i]=0;
   vx[i]=0;
   vy[i]=0;
   energ[i]=0;
   err[i]=0;
}

x[0]= 1-e;
y[0]=0;
vx[0]=0;
vy[0]=sqrt((1.0 +e)/(1.0 -e));

energ[0]=Energia(x[0],y[0],vx[0],vy[0]);

for(int i=0; i<pasos-1; i++){  //x explícita e y implícita

x[i+1]=x[i]+h*vx[i]; //exp
y[i+1]=y[i]+h*vy[i]; //imp

vx[i+1]=vx[i]+h*f_Kepler(x[i+1],y[i+1]);
vy[i+1]=vy[i]+h*f_Kepler(y[i+1],x[i+1]);




energ[i]=Energia(x[i],y[i],vx[i],vy[i]);
err[i]=Calcerror(x[i],y[i]);

//printf("Paso %d. X=%f  Y= %f Energia=%f  INI X0 = %f Y0 =%f\n", i, x[i], y[i], energ[i], x[0],y[0]);

}
printf("Valores iniciales X = %f Y = %f", x[0],y[0]);
printEsymp(x,y,energ,err,pasos);

}

void printEexplicito(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f1;
 f1=fopen("Euler_explicito00001.txt", "w");
	if(f1==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.0001;

            fprintf(f1,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f1);



}

void printEsymp(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f2;
 f2=fopen("Euler_symp.txt", "w");
	if(f2==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.001;
            fprintf(f2,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);

    }
    fclose(f2);


}


void printMidpoint(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f3;
 f3=fopen("Midpoint1.txt", "w");
	if(f3==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
             aux=i*0.05;
            fprintf(f3,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f3);

}


void printVerlet(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f4;
 f4=fopen("Verlet.txt", "w");
	if(f4==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
             aux=i*0.05;
            fprintf(f4,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f4);
}


float Calcerror(float x, float y){
float phi,r,r_teo,err,d;
float pi;
pi=4*atan(1.0);
d=1-e*e;
if(x>0 && y>=0) phi=(float)atan(y/x);
if(x==0 && y>0) phi=0.5*pi;
if(x>0 && y <0) phi=(float)atan(y/x) + 2*pi;
if(x<0) phi=(float)atan(y/x) + pi;
if(x==0 && y<0) phi=1.5*pi;

r=sqrt(x*x+y*y);
r_teo= (d)/(1.0+e*(float)cos(phi));

err=(r_teo-r)*(r_teo - r);
return err;

}

void CalculoAnalitico(int pasos){
    float pi;
pi=4*atan(1.0);
float delta=(2*pi)/(float)pasos;
//printf("delta %f \n", delta);
float d=1-e*e;
float x[pasos];
float y[pasos];

float phi,r;
phi=0;

for(int i=0; i<pasos; i++){
r= (d)/(1.0+e*(float)cos(phi));
x[i]=r*cos(phi);
y[i]=r*sin(phi);
phi=phi+delta;
//printf("Paso %d. phi = %f, r = %f,  x = %f , y= %f \n", i,phi, r, x[i],y[i]);
//getchar();
}

FILE *f5;
 f5=fopen("Analitico.txt", "w");
	if(f5==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){

            fprintf(f5,"%f   %f\n",x[i],y[i]);
    }
    fclose(f5);

}


