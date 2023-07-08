# include <iostream>
# include <cmath>
# include <fstream>
# include "gsl_rng.h" //Libreria para generación de números aleatorios

using namespace std;

gsl_rng *tau;           // Inicializo el puntero para generar números aleatorios

int main()
{
    int i = 0;              // Variable contador
    int j = 0;              // Variable contador
    int k = 0;              // Variable contador
    int l = 0;              // Variable contador
    int t = 0;              // Variable contador general
    int w = 0;              // Variable contador de patrones
    int N = 20;            // Variable para determinar la dimensión de la matriz neuronal
    int fila = 0;           // Variable para determinar la posición de la neurona estudiada
    int columna = 0;        // Variable para determinar la posición de la neurona estudiada
    int patron = 50;       // Variable para indicar el número de patrones máximo a incluir en el código
    int ncuadrado = N * N;  // Variable para evitar al cálculo de N*N en cada paso
    int neuronas[N][N];            // Declaro la matriz de neuronas
    int patrones[N][N][patron];  // Cubo para guardar las diversas configuraciones 
    double auxiliarpat[N][N][patron];   // Cubo para guardar las configuraciones modificadas
    double tetagm[N][N];       // Matriz para guardar los tetas de todos los patrones
    double T = 0.0001;         // Temperatura del sistema
    double proppatron[patron];      // Proporción de neuronas encendidas en cada patron
    double energia = 0.0;           // Variable para almacenar el cambio de energia entre estados del sistema
    double epsilon = 0.0;           // Variable para guardar un número aleatorio
    double solapamiento[patron];    // Para calcular el solapamiento de cada patrón con la matriz en cada paso montecarlo
    extern gsl_rng *tau;    //Puntero al estado del número aleatorio
    int semilla = 135254;     //Semilla del generador de números aleatorios
    tau = gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    ofstream solapamientos;    // Declaramos un fichero para escribir la matriz de neuronas

    solapamientos.open("solappokemongeneral.txt");

    for (t = 0; t <= (patron - 1); t++)
    {
        // Genero las matrices de neuronas que luego quiero recordar

        for (w = 0; w <= t; w++)
        {
            for (i = 0; i <= (N - 1); i++)
            {
                for (j = 0; j <= (N - 1); j++)
                {
                    patrones[i][j][w] = gsl_rng_uniform_int(tau,2);  //número aleatorio entero [0,1]
                }
            }
        }
        
        // Genero las matrices de neuronas inicial         
    
        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                neuronas[i][j] = gsl_rng_uniform_int(tau,2);  //número aleatorio entero [0,1]
            }
        }

        // Calculamos el valor de a^mu para todos los patrones

        for (w = 0; w <= t; w++)
        {
            proppatron[w] = 0.0;

            for (i = 0; i <= (N - 1); i++)
            {
                for (j = 0; j <= (N - 1); j++)
                {
                    proppatron[w] = proppatron[w] + patrones[i][j][w];
                } 
            }

            proppatron[w] = proppatron[w]/(1.0 * ncuadrado); 
        }

        // Voy a crear una matriz auxiliar en la que guardo los valores (epsilon - a^mu) para cada patrón
    
        for(w = 0; w <= t; w++)
        {
            for (i = 0; i <= (N - 1); i++)
            {
                for (j = 0; j <= (N - 1); j++)
                {
                    auxiliarpat[i][j][w] = patrones[i][j][w] - proppatron[w];
                
                    if (w == 0)
                    {
                        tetagm[i][j] = 0.0;
                    }
                }
            }
        }
        // Calculo la matriz de tetas correspondiente a la unión de todos los patrones

        for (i = 0; i <= (N - 1); i++)
        {
            for (j = 0; j <= (N - 1); j++)
            {
                for (w = 0; w <= t; w++)
                {
                    for (k = 0; k <= (N - 1); k++)
                    {
                        for (l = 0; l<= (N - 1); l++)
                        {
                            if ((i != k) || (j != l))
                            {
                                tetagm[i][j] = tetagm[i][j] + auxiliarpat[i][j][w] * auxiliarpat[k][l][w];
                            }
                        }
                    }
                }

                tetagm[i][j] = tetagm[i][j]/(2.0);      // No divido por N**2 porque lo hago luego
            }
        } 

        // Inicializo el bucle de iteraciones para poner en práctica el algoritmo de metrópolis

        for (i = 0; i <= 50 * ncuadrado; i++)
        {
            if (i == 50 * ncuadrado)
            {
                for (w = 0; w <= t; w++)
                {
                    solapamiento[w] = 0.0;

                    for (k = 0; k <= (N - 1); k++)
                    {
                        for (l = 0; l <= (N - 1); l++)
                        {
                            solapamiento[w] = solapamiento[w] + auxiliarpat[k][l][w] * (neuronas[k][l] - proppatron[w]);
                        }
                    }

                    solapamiento[w] = solapamiento[w]/(ncuadrado * proppatron[w] * (1 - proppatron[w]));
                    solapamientos << solapamiento[w] << endl;

                    if (w == t)
                    {
                        solapamientos << endl;
                    }
                }
            }

            // Escojo un punto al azar del sistema

            fila = gsl_rng_uniform_int(tau, N);  //número aleatorio entero [0,N - 1]
            columna = gsl_rng_uniform_int(tau, N);  //número aleatorio entero [0, N - 1]

            // Reinicio el valor de la variable de energía (incremento de energía)

            energia = 0.0;

            // Calculo el incremento de energía

            for (w = 0; w <= t; w++)
            {
                for (k = 0; k <= (N - 1); k++)
                {
                    for (l = 0; l <= (N - 1); l++)
                    {
                        if ((fila != k) || (columna != l))
                        {
                            energia = energia + auxiliarpat[fila][columna][w] * auxiliarpat[k][l][w] * neuronas[k][l];
                        }
                    }
                }
            }

            energia = (energia - tetagm[fila][columna])/(ncuadrado * 1.0);

            // Determino si el paso se realiza de 1 a 0 o viceversa

            if (neuronas[fila][columna] == 0)
            {
                energia = energia * (-1.0);
            }

            // Calculo el valor de p

            energia = exp(energia/(-1.0 * T));

            if (energia > 1)
            {
                energia = 1;
            }

            // Generamos un número aleatorio entre 0 y 1 y comparamos con energia. Si es menor cambiamos el estado de la neurona

            epsilon = gsl_rng_uniform(tau);  //número aleatorio real en [0,1]

            if (epsilon < energia)
            {
                neuronas[fila][columna] = 1 - neuronas[fila][columna];
            }
        }
    }    

    solapamientos.close();
}