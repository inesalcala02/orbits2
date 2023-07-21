#include <iostream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <vector>

const double G = 0.000295912208286; ///DISTANCIA EN u.a. (1 u.a. = 149 597 870 km)
const double m0 = 1.00000597682; ///TOMASMOS MASA DEL SOL = 1, M0 TIENE EN CUENTA
                                ///LOS PLANETAS PRÓXIMOS AL SOL



struct Body{

double x;
       double y;
       double z;

       double vx;
       double vy;
       double vz;

       double mass;



};

void InitializeSolarSystem();
double KeplerForce_x( double x, double x1, double x2, std::vector <Body> bodies , Body planet);
double KeplerForce_y( double x, double x1, double x2, std::vector <Body> bodies , Body planet);
double KeplerForce_z( double x, double x1, double x2, std::vector <Body> bodies , Body planet);
double energia( std::vector <double> x,std::vector <double> y,std::vector <double> z,std::vector <double> vx,std::vector <double> vy,std::vector <double> vz );

///metodos
void Explicit(int steps, double h);
void Verlet(int steps, double h);


std::vector <Body> planets;

bool doVerlet = true;
bool doExplicit = false;
bool doSymplectic = false;
bool doRK4 = false;

int main(){

///Inicializar
InitializeSolarSystem();

 ///VERLET (sale)
    if( doVerlet )
    {
        Verlet(1000,200);
        puts("Metodo de Verlet acabado ... ");
        getchar();
    }

    ///EXP EULER
    if( doExplicit )
    {
        Explicit(20000, 10);
        puts("Metodo exp acabado ... ");
        getchar();
    }

    ///SYMPLECTIC EULER (sale)
    if( doSymplectic )
    {
       // SympcEuler(2000,100);
        puts("Metodo symplectic Euler acabado ... ");
        getchar();
    }

    ///RK4
    if( doRK4 )
    {
       // RK4(1000,200);
        puts("Metodo RK4 acabado ... ");
        getchar();
    }




return 0;
}

void Explicit(int steps, double h)
{
FILE *f1 = fopen("Exp_Jup.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f2 = fopen("Exp_Sat.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f3 = fopen("Exp_Ur.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f4 = fopen("Exp_Nep.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f5 = fopen("Exp_Pl.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

std::vector <double> x, y, z, vx, vy, vz;
double energy;

int aux;
int i=0;

   while (i<steps)
    {  //printf("Paso nuevo %d, h vale %lf\n", i, h);
        for ( auto &b : planets )
        {
          x.push_back( b.x + b.vx * h);                              //para guardar datos
          vx.push_back( b.vx + h * KeplerForce_x(b.x,b.y,b.z, planets, b));
       //   printf("Kepler force en x = %lf\n",KeplerForce_x(b.x, planets, b) );
          y.push_back( b.y + b.vy * h);
          vy.push_back( b.vy + h * KeplerForce_y(b.y,b.x,b.z, planets, b));

          z.push_back( b.z + b.vz * h);
          vz.push_back( b.vz + h * KeplerForce_z(b.z,b.x,b.y, planets, b));

        //if(b.mass == 0.000954786104043)  printf("Planeta nuevo: x = %f , vx = %f , y = %f , vy = %f , z = %f, vz = %f\n", b.x, b.vx, b.y,b.vy,b.z,b.vz);
          //getchar();
        }
    aux= i*planets.size();
    for(int j=0; j<planets.size(); j++)
    {
        planets[j].x = x[j+aux];
        planets[j].vx = vx[j+aux];
        planets[j].y = y[j+aux];
        planets[j].vy = vy[j+aux];
        planets[j].z = z[j+aux];
        planets[j].vz = vz[j+aux];
    }

   // energy= energia(x,y,z,vx,vy,vz);

    fprintf(f1,"%lf   %lf  %lf   %lf\n", h * (double) i, x[aux], y[aux], z[aux]);
    fprintf(f2,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+1], y[aux+1], z[aux+1]);
    fprintf(f3,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+2], y[aux+2], z[aux+2]);
    fprintf(f4,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+3], y[aux+3], z[aux+3]);
    fprintf(f5,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+4], y[aux+4], z[aux+4]);


    i++;

    }
}

void Verlet(int steps, double h)
{
 FILE *f1 = fopen("Ver_Jup.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f2 = fopen("Ver_Sat.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f3 = fopen("Ver_Ur.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f4 = fopen("Ver_Nep.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

FILE *f5 = fopen("Ver_Pl.txt", "w");
if(!f1) { printf("Error al abrir el archivo de texto Exp_sist."); exit(1); }

std::vector <double> x, y, z, vx, vy, vz,vx_m, vy_m, vz_m;


int aux;
int i=0;

   while (i<steps)
    {  //printf("Paso nuevo %d, h vale %lf\n", i, h);
        for ( auto &b : planets )
        {
            vx_m.push_back(b.vx + 0.5 * h * KeplerForce_x(b.x,b.y,b.z, planets, b));
            x.push_back( b.x + vx_m.back() * h);

            vy_m.push_back ( b.vy + 0.5 * h * KeplerForce_y(b.y,b.x,b.z, planets, b));
            y.push_back( b.y + vy_m.back() * h);

            vz_m.push_back( b.vz + 0.5 * h * KeplerForce_z(b.z,b.x,b.y, planets, b));
            z.push_back( b.z + vz_m.back() * h);

        }

        aux= i*planets.size();
    for(int j=0; j<planets.size(); j++)
    {
        planets[j].x = x[j+aux];
        planets[j].y = y[j+aux];
        planets[j].z = z[j+aux];
    }
        int a=0;
        for ( auto &b : planets ) ///PARA ACTUALIZAR LAS VELOCIDADES CON TODAS LAS POSICIONES NUEVAS
        {
           vx.push_back( vx_m[aux+a] + h * 0.5 * KeplerForce_x(b.x,b.y,b.z, planets, b));
           vy.push_back( vy_m[aux+a] + h * 0.5 * KeplerForce_y(b.y,b.x,b.z, planets, b));
           vz.push_back( vz_m[aux+a] + h * 0.5 * KeplerForce_z(b.z,b.x,b.y, planets, b));
           a++;
        }

        for(int j=0; j<planets.size(); j++)
    {
        planets[j].vx = vx[j+aux];
        planets[j].vy = vy[j+aux];
        planets[j].vz = vz[j+aux];
    }




   // energy= energia(x,y,z,vx,vy,vz);

    fprintf(f1,"%lf   %lf  %lf   %lf\n", h * (double) i, x[aux], y[aux], z[aux]);
    fprintf(f2,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+1], y[aux+1], z[aux+1]);
    fprintf(f3,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+2], y[aux+2], z[aux+2]);
    fprintf(f4,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+3], y[aux+3], z[aux+3]);
    fprintf(f5,"%lf   %lf  %lf   %lf\n", h * (float) i, x[aux+4], y[aux+4], z[aux+4]);


    i++;

    }
}




double KeplerForce_x( double x, double x1, double x2, std::vector <Body> bodies , Body planet)
{

       double force_x ;
       double m = planet.mass;
       double dist = pow(x*x + x1*x1 + x2*x2,1.5);
        force_x = -( m0 ) * x / (double)dist; ///distancia al origen (sol)


 //if(planet.mass == 0.000954786104043)printf("ESTAMOS EN LA FUNCION DE KEPLERFORCE\n Tenemos x = %lf masa = %lf G = %lf fuerza_0 = %lf  m0 = %lf\n", x, m,G, force_x, m0);

        for ( auto &b : bodies )   ///como no utilizar el planeta para el cual queremos calcular la fuerza??
        {
               double x_p = b.x;
               double y_p = b.y;
               double z_p = b.z;
               double m_p = b.mass;

              if(m_p != m)
                {   double cubo = pow((x-x_p)*(x-x_p) + (x1-y_p)*(x1-y_p) + (x2-z_p)*(x2-z_p),1.5);
                    force_x +=  m_p * (x_p - x)/(double)cubo;
   //           if(planet.mass == 0.000954786104043)  printf("m_p = %lf y m = %lf entonce suma %lf  (x_p =%lf y x =%lf)\n", m_p, m, m_p * (x_p - x)/(double)cubo,x_p,x);
                }
        }
     //if(planet.mass == 0.000954786104043)   printf("DESPUES DE LA SUMA DE TODOS LOS PLANETAS LA FUERZA ES %lf\n",force_x );
     force_x = G * force_x;
   //if(planet.mass == 0.000954786104043) printf("DESPUES DE MULTIPLICAR G LA FUERZA ES %lf\n",force_x );
     return force_x;
}

double KeplerForce_y( double x, double x1, double x2, std::vector <Body> bodies , Body planet)
{

       double force_x ;
       double m = planet.mass;
       double dist = pow(x*x + x1*x1 + x2*x2,1.5);
        force_x = -( m0 ) * x / (double)dist;
         ///distancia al origen (sol)


 //printf("ESTAMOS EN LA FUNCION DE KEPLERFORCE\n Tenemos x = %lf masa = %lf G = %lf fuerza_0 = %lf  m0 = %lf\n", x, m,G, force_x, m0);

        for ( auto &b : bodies )   ///como no utilizar el planeta para el cual queremos calcular la fuerza??
        {
               double x_p = b.x;
               double y_p = b.y;
               double z_p = b.z;
               double m_p = b.mass;

              if(m_p != m)
                {   double cubo =  pow((x-y_p)*(x-y_p) + (x1-x_p)*(x1-x_p) + (x2-z_p)*(x2-z_p),1.5);
                    force_x +=  m_p * (x_p - x)/(double)cubo;
              //  printf("m_p = %lf y m = %lf entonce suma %lf  (x_p =%lf y x =%lf)\n", m_p, m,force_x,x_p,x);
                }
        }
     force_x = G * force_x;
  //  printf("DESPUES DE LA SUMA DE TODOS LOS PLANETAS LA FUERZA ES %lf\n",force_x );
     return force_x;
}

double KeplerForce_z( double x, double x1, double x2, std::vector <Body> bodies , Body planet)
{

       double force_x ;
       double m = planet.mass;
        double dist = pow(x*x + x1*x1 + x2*x2,1.5);
        force_x = -( m0 ) * x / (double)dist;
         ///distancia al origen (sol)


 //printf("ESTAMOS EN LA FUNCION DE KEPLERFORCE\n Tenemos x = %lf masa = %lf G = %lf fuerza_0 = %lf  m0 = %lf\n", x, m,G, force_x, m0);

        for ( auto &b : bodies )   ///como no utilizar el planeta para el cual queremos calcular la fuerza??
        {
                double x_p = b.x;
               double y_p = b.y;
               double z_p = b.z;
               double m_p = b.mass;

              if(m_p != m)
                {   double cubo =  pow((x-z_p)*(x-z_p) + (x1-x_p)*(x1-x_p) + (x2-y_p)*(x2-y_p),1.5);
                    force_x +=  m_p * (x_p - x)/(double)cubo;
               // printf("m_p = %lf y m = %lf entonce suma %lf  (x_p =%lf y x =%lf)\n", m_p, m,force_x,x_p,x);
                }
        }
     force_x = G * force_x;
   // printf("DESPUES DE LA SUMA DE TODOS LOS PLANETAS LA FUERZA ES %lf\n",force_x );
     return force_x;
}


void InitializeSolarSystem()
{

// Jupiter
Body Jupyter;
Jupyter.x = -3.5023653;
Jupyter.y = -3.8169847;
Jupyter.z = -1.5507963;
Jupyter.vx = 0.00565429;
Jupyter.vy = -0.00412490;
Jupyter.vz = -0.00190589;
Jupyter.mass = 0.000954786104043;

planets.push_back(Jupyter);

Body Saturn;
Saturn.x = 9.0755314;
Saturn.y = -3.0458353;
Saturn.z = -1.6483708;
Saturn.vx = 0.00168318;
Saturn.vy = 0.00483525;
Saturn.vz = 0.00192462;
Saturn.mass = 0.000285583733151;

planets.push_back(Saturn);

Body Uranus;
Uranus.x = 8.3101420;
Uranus.y = -16.2901086;
Uranus.z = -7.2521278;
Uranus.vx = 0.00354178;
Uranus.vy = 0.00137102;
Uranus.vz = 0.00055029;
Uranus.mass = 0.0000437273164546;

planets.push_back(Uranus);

Body Neptune;
Neptune.x = 11.4707666;
Neptune.y = -25.7294829;
Neptune.z = -10.8169456;
Neptune.vx = 0.00288930;
Neptune.vy = 0.00114527;
Neptune.vz = 0.00039677;
Neptune.mass = 0.0000517759138449;

planets.push_back(Neptune);

Body Pluto;
Pluto.x = -15.5387357;
Pluto.y = -25.2225594;
Pluto.z = -3.1902382;
Pluto.vx = 0.00276725;
Pluto.vy = -0.00170702;
Pluto.vz = -0.00136504;
Pluto.mass = 0.00000000769230769;

planets.push_back(Pluto);


}
