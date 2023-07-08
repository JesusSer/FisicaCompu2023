# include <stdio.h>
# include <math.h>
# include "gsl_rng.h" //Libreria para generación de números aleatorios

gsl_rng *tau;

int main()
{
    int i,j,k;  //Variables auxiliares
    double x,y,z;  //Variables auxiliares
    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=135254; //Semilla del generador de números aleatorios
    
    
    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla
    
    for (i=1;i<=10;i++)
    {
        x=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
        j=gsl_rng_uniform_int(tau,20);  //número aleatorio entero [0,19]
        printf("%lf\t%i\n",x,j);
    }
}