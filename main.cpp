#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <vector>


const double H0 = -0.5; // Energía inicial, en teoría debería conservarse
const double phi1 = 0;  // ángulo phi*
const double e = 0.6; //excentricidad

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
void printMidpoint(std::vector <float> x, std::vector<float> y, std::vector<float> en, std::vector<float> err);
void printVerlet(float x[], float y[],float ener[],float err[], int pasos);
float Calcerror(float x, float y); ///El valor de d hay que ponerlo a mano
void CalculoAnalitico(int pasos);

void RK4(float pasos, float h);
// Es mejor omitir, para tener claro los miembros que vienen de std::
//using namespace std;

bool doVerlet = true;
bool doMidPoint = false;
bool doSymplectic = false;
bool doRK4 =true;
int main()
{

    //definición de d y L0
    float d = 1.0 - e*e;
    float L0 = sqrt(d);

    printf("Parametros del sistema: H0 = %f  e = %f L0 = %f d = %f\n",H0,e, L0, d);
    puts("Pulsa Intro para continuar ... ");
    getchar();

    // crear array de los pasos temporales para distintos métodos

    float h[4];
    int steps[4];
    h[0]= 0.0005;
    steps[0]=400000 ;//euler explicito
    h[1]=0.01;
    steps[1]=40000; //symplectic Euler, implicit midpoint, Verlet


    ///SOLUCION ANALITICA (sale)
    //CalculoAnalitico(steps[1]);
    puts("Calculo analitico acabado ... ");

    ///EULER EXPLÍCITO (sale)
    //EulerExplicito(steps[0], h[0]);
    // puts("Metodo Euler explicito acabado ... ");
    //getchar();

    ///VERLET (sale)
    if( doVerlet )
    {
        Verlet(steps[1],h[1]);
        puts("Metodo de Verlet acabado ... ");
        getchar();
    }

    ///IMPLICIT MIDPOINT
    if( doMidPoint )
    {
        ImpMidpoint(steps[1], h[1]);
        puts("Metodo implicit midpoint acabado ... ");
        getchar();
    }

    ///SYMPLECTIC EULER (sale)
    if( doSymplectic )
    {
        SympcEuler(steps[1],h[1]);
        puts("Metodo symplectic Euler acabado ... ");
        getchar();
    }

    ///RK4
    if( doRK4 )
    {
    RK4(steps[1],h[1]);
        puts("Metodo RK4 acabado ... ");
        getchar();
    }


    return 0;
}


void EulerExplicito(int pasos, float h){


// ---> No es necesario definir el numero de pasos cuando se usa std::vector
float x[pasos];
float y[pasos];
float vx[pasos];
float vy[pasos];
float energ[pasos];
float err[pasos];

// ---> tampoco inicializar, puesto que el vector estara vacio
for (int i=0; i<pasos; i++){  //inicializar vectores
   x[i]=0;
   y[i]=0;
   vx[i]=0;
   vy[i]=0;
   energ[i]=0;
   err[i]=0;
}

// Habria que usar esto, y push_back
// std::vector <float> x,y,vx,vy,energ,err;

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
err[i]=fabs(-0.5 - energ[i])/(float)energ[i];

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
err[i+1]=(-0.5 - energ[i])/(-0.5);
//printf("Paso %d \n", i);

}

printVerlet(x,y,energ,err,pasos);
}



void ImpMidpoint(int pasos, float h){

    // Setting vectors and initial values
    std::vector <float> x, y, vx, vy;
    std::vector <float> energy, err;

    x.push_back(1-e);
    y.push_back(0.01); // Por que 0.01 y no 0?
    vx.push_back(0);
    vy.push_back(sqrt((1.0 +e)/(1.0 -e)));
    energy.push_back(-0.5);
    err.push_back(0);

    for(int i=0; i<pasos-1;i++){

        /// Getting the last vector element
        float Xo = x.back();
        float Yo = y.back();

        float VXo = vx.back();
        float VYo = vy.back();

        float k1_vx = f_Kepler(Xo, Yo);
        float k1_vy = f_Kepler(Yo, Xo);

        float k1_x = Xo + h * (VXo + 0.5 * h * k1_vx);
        float k1_y = Yo + h * (VYo + 0.5 * h * k1_vy);

        float k2_vx = f_Kepler( k1_x, k1_y );
        float k2_vy = f_Kepler( k1_y, k1_x );

        x.push_back( Xo + h * VXo + 0.5 * h * h * k1_vx );
        y.push_back( Yo + h * VYo + 0.5 * h * h * k1_vy );

        vx.push_back( VXo + 0.5 * h * (k1_vx + k2_vx) );
        vy.push_back( VYo + 0.5 * h * (k1_vy + k2_vy) );

        energy.push_back( Energia( x.back(), y.back(), vx.back(), vy.back() ) );
        err.push_back( Calcerror( x.back(), y.back() ) );

        //printf("Paso %d. X=%f  Y= %f Energia=%f  INI X0 = %f Y0 =%f\n", i, x[i], y[i], energ[i], x[0],y[0]);
        //getchar();
    }

    FILE *f3 = fopen("Midpoint1.txt", "w");
    if(!f3) { printf("Error al abrir el archivo de texto."); exit(1); }

    for(unsigned int i=0; i < x.size(); i++)
        fprintf(f3,"%f   %f    %f   %f   %f\n", 0.05 * (float) i, x[i], y[i], energy[i], err[i]);
    fclose(f3);
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
err[i]=-(-0.5 - energ[i])/(-0.5);

//printf("Paso %d. X=%f  Y= %f Energia=%f  INI X0 = %f Y0 =%f\n", i, x[i], y[i], energ[i], x[0],y[0]);

}
printf("Valores iniciales X = %f Y = %f", x[0],y[0]);
printEsymp(x,y,energ,err,pasos);

}



void RK4(float pasos, float h){

std::vector <float> x, y, vx, vy;
    std::vector <float> energy, err;

    x.push_back(0.4);
    y.push_back(0.0);
    vx.push_back(0);
    vy.push_back(2);
    energy.push_back(-0.5);
    err.push_back(0.0);

 for(int i=0; i<pasos-1;i++){

        /// Getting the last vector element
        float Xo = x.back();
        float Yo = y.back();

        float VXo = vx.back();
        float VYo = vy.back();

        //printf( "xo = %f,  vx0 = %f, y0 = %f, vy0 = %f \n", Xo, VXo, Yo, VYo);

        float k1_vx = f_Kepler(Xo, Yo);
        float k1_vy = f_Kepler(Yo, Xo);

        float k1_x = VXo;
        float k1_y = VYo;


        //printf( "VX k1 = %f, VY k1 = %f \n", k1_vx,k1_vy);


        float k2_vx = f_Kepler(Xo + 0.5 * h * k1_x, Yo + 0.5 * h * k1_y );
        float k2_vy = f_Kepler(Yo + 0.5 * h * k1_y, Xo + 0.5 * h * k1_x );

        float k2_x = VXo + h * 0.5 * k1_vx;
        float k2_y = VYo + h * 0.5 * k1_vy;


          //printf( "VX k2 = %f, VY k2 = %f \n", k2_vx,k2_vy);

        float k3_vx = f_Kepler(Xo + 0.5 * h * k2_x, Yo + 0.5 * h * k2_y );
        float k3_vy = f_Kepler(Yo + 0.5 * h * k2_y, Xo + 0.5 * h * k2_x );

        float k3_x = VXo + h * 0.5 * k2_vx;
        float k3_y = VYo + h * 0.5 * k2_vy;

         // printf( "VX k3 = %f, VY k3 = %f \n", k3_vx,k3_vy);

        float k4_vx = f_Kepler(Xo + h * k3_x, Yo + h * k3_y );
        float k4_vy = f_Kepler(Yo + h * k3_y, Xo + h * k3_x );

        float k4_x = VXo + h * k3_vx;
        float k4_y = VYo + h * k3_vy;

          //printf( "VX k4 = %f, VY k4 = %f \n", k4_vx,k4_vy);


        vx.push_back( VXo + h * (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx) * 0.1666666 );
        vy.push_back( VYo + h * (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy) * 0.1666666 );


        x.push_back( Xo +  h * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x) * 0.1666666 );
        y.push_back( Yo +  h * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y) * 0.1666666 );



        energy.push_back( Energia( x.back(), y.back(), vx.back(), vy.back() ) );
        err.push_back(-(-0.5 - energy[i])/(-0.5));

        //printf("Paso %d. X=%f  Y= %f Energia=%f  VEL VX = %f VY =%f\n", i, x[i], y[i], energy[i], vx[i],vy[i]);
       // getchar();
    }

    printf("X=%f Y =%f VX=%f VY=%f", x[pasos-1], y[pasos-1], vx[pasos-1],vy[pasos-1]);
    FILE *f5 = fopen("RK4_prueba.txt", "w");
    if(!f5) { printf("Error al abrir el archivo de texto."); exit(1); }

    for(unsigned int i=0; i < x.size(); i++)
        fprintf(f5,"%f   %f    %f   %f   %f\n", 0.01 * (float) i, x[i], y[i], energy[i], err[i]);
    fclose(f5);
}







void printEexplicito(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f1;
 f1=fopen("Euler_explicito.txt", "w");
	if(f1==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.05;

            fprintf(f1,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);
    }
    fclose(f1);



}

void printEsymp(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f2;
 f2=fopen("Euler_symp_0001.txt", "w");
	if(f2==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
            aux=i*0.001;
            fprintf(f2,"%f   %f    %f   %f   %f\n",aux,x[i],y[i],ener[i], err[i]);

    }
    fclose(f2);


}




void printVerlet(float x[], float y[],float ener[],float err[], int pasos){
    float aux;
FILE *f4;
 f4=fopen("Verlet_prueba.txt", "w");
	if(f4==NULL){printf("Error al abrir el archivo de texto.");exit(1);	}

    for(int i=0;i<pasos;i++){
             aux=i*0.01;
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


