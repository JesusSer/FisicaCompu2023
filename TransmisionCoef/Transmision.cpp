// Aquí vamos a desarrollar el código para resolver la ecuación de Schrödinger 

// Discretizamos espacio y tiempo: x_j = jh y t_n = ns, con j yendo de 0 a N
// La función de onda en cada punto del retículo creado es phi_j,n 
// Condiciones de contorno: phi_0,n = phi_N,n = 0. Esto implica X_N = X_0 = 0

// Usamos el algoritmo phi_j,n+1 = frac(1-(isH_d)/2)(1+(isH_d)/2) * phi_j,n
// Llegamos a un sistema del tipo X_j+1,n + [-2+2i/s-V_jgorrito]X_j,n + X_j-1,n = 4*i*phi_j,n/sgorrito
// Tomamos los cambios sgorrito = s/h^2 y V_jgorrito = h^2*V_j

// Aparecerán en la resolución del sistema triangular otras constantes alpha y beta
// Alpha es constante en el tiempo
// Tomamos como función de onda inicial una onda con amplitud gaussiana
// La ecuación de la gaussiana es phi(x,0) = e^(i(k_0)x)e^(-(x-x_0)^2/(2*sigma^2))

// Debemos dar un parámetro nciclos para obtener k_0 de modo que (k_0)*N*h = 2*pi*nciclos
// nciclos va de 1 a N/4, así un ciclo tendrá 4 puntos como mínimo
// Podemos dar como valores iniciales: x_0 = N*h/4 y sigma = H*h/16, aunque es recomendable escribir la función de onda
// inicial de maner genérica y jugar con x_0 y sigma

// El potencial que usaremos tendrá una anchura de N/5 y estará centrado en N/2, con una altura proporcional ala energía de la onda
// incidente, por ejemplo, altura = lambda * (k_0)^2, con lambda = 0.3

// Así, la función de onda en el retículo quedará phi(x,0) = e^(i(k_0gorrito)j)e^(-8(4j-N)^2/N^2), con k_0gorrito = (k_0)*h
// El potencial queda V_jgorrito = (V_j)*h= [0] si j no está en (2N/5, 3N/5) ó lambda * (k_0gorrito)^2 en el resto de los casos

// Tomaremos sgorrito = 1/(4*(k_0gorrito)^2)

// Así, al principio solo tenemos que definir N, nciclos y lambda

// El algoritmo para generar la onda será:
// 1 -> Doy los parámetros N, nciclos y lambda. Con ellos genero sgorrito, k_0gorrito, V_jgorrito, phi_j,0, condiciones de contorno y alpha
// 2 -> Calculo beta
// 3 -> Calculo X
// 4 -> Calculo phi_j,n+1
// 5 -> n = n+1, vuelvo a 2

// -------------------- EMPIEZO CON EL CÓDIGO ---------------

// Incluyo las librerías necesarias

# include <iostream>
# include <cmath>
# include <fstream>
# include "complex.h" //Biblioteca externa para trabajar con numeros complejos
# include "gsl_rng.h" //Libreria para generación de números aleatorios

using namespace std;

gsl_rng *tau;           // Inicializo el puntero para generar números aleatorios

int main()
{
    int i = 0;                  // Variable contador
    int t = 0;                  // Variable contador del tiempo
    int N = 2000;                // Número de divisiones del espacio
    int tiempomaximo = 0;       // Tiempo en que se alcanza el primer máximo relativo en el detector derecho
    int particulas = 0;         // Número de partículas que atraviesan la barrera de potencial

    double valorlambda = 10.0;   // Valor de lambda
    double vgorrito[N + 1];     // Vector para guardar la barrera de potencial
    double norma = 0.0;         // Variable para almacenar la norma de la función de onda y normalizar
    double nciclos = N / 4.0;   // Variable para almacenar nciclos  
    double vectoronda[N + 1];   // Vector para almacenar la norma al cuadrado de cada posición de la onda
    double kgorrito;            // Variable para guardar el valor de k
    double sgorrito;            // Variable para guardar el valor de s
    double valorvgorrito;       // Variable para guardar el valor constante del pozo de potencial
    double liminf;              // Variable auxiliar para rellenar el vector de potencial
    double limsup;              // Variable auxiliar para rellenar el vector de potencial
    double probinteres = 0.0;   // Probabilidad de encontrar la partícula en el detector derecho
    double variableauxiliar = 0.0;  // Variable auxiliar para realizar diversos cálculos
    double aleatorio;           // Variable que almacena el número aleatorio generado en cada caso
    double transmision;         // Variable para guardar el coeficiente de transmisión obtenido

    bool condicion = true;      // Condición para determinar el máximo relativo de la probabilidad en el detector a la derecha

    extern gsl_rng *tau;        // Puntero al estado del número aleatorio
    int semilla = 1345249;       // Semilla del generador de números aleatorios
    tau = gsl_rng_alloc(gsl_rng_taus); // Inicializamos el puntero
    gsl_rng_set(tau,semilla);   // Inicializamos la semilla

    fcomplex phi[N + 1];        // Vector para almacenar el valor de la función de onda en cada posición
    fcomplex vectoralpha[N];    // Vector alpha auxiliar para calcular la función de onda en cada paso
    fcomplex vectorbeta[N];     // Vector beta auxiliar para calcular la función de onda en cada paso
    fcomplex vectorchi[N + 1];  // Vector chi auxiliar para calcular la función de onda en cada paso

    // Calculo el resto de parámetros

    kgorrito = 2 * (M_PI) * nciclos / (N * 1.0);
    sgorrito = 1.0 / (4.0 * kgorrito * kgorrito);
    valorvgorrito = valorlambda * kgorrito * kgorrito;

    // Relleno el vector del potencial

    liminf = 2.0 * N / 5.0;
    limsup = 3.0 * N / 5.0;

    for (i = 0; i <= N; i++)
    {
        if ((liminf <= i) && (i <= limsup))
        {
            vgorrito[i] = valorvgorrito;
        }
        else
        {
            vgorrito[i] = 0.0;
        }
    }

    // Calculo la función de onda inicial. Comenzamos aplicando las condiciones contorno

    phi[0] = Complex(0.0, 0.0);
    phi[N] = Complex(0.0, 0.0);

    for (i = 1; i < N; i++)
    {
        phi[i] = (RCmul(exp((-(8.0) * (4 * i - N) * (4 * i - N) / (N * N))), (Complex(cos(kgorrito * i), sin(kgorrito * i)))));
        vectoronda[i] = Cabs(phi[i]) * Cabs(phi[i]);
        norma = norma + vectoronda[i];
    }

    for (i = 1; i < N; i++)
    {
       phi[i] = RCmul(1.0 / sqrt(norma) , phi[i]); 
    }

    // Relleno ahora el valor del vector auxiliar alpha. Empiezo por el final

    vectoralpha[N - 1] = Complex(0.0, 0.0);

    for (i = 2; i <= N; i++)
    {
        vectoralpha[N - i] = Cdiv(Complex(-1.0, 0.0), Cadd(Complex(-2.0 - vgorrito[N - i + 1], 2.0 / sgorrito), vectoralpha[N - i + 1]));
    }

    // Finalmente aplicamos el algoritmo para un tiempo discretizado
    
    for (t = 0; t <= 5000; t++) 
    {
        // Calculamos el vector de betas. Empezamos por el final

        vectorbeta[N - 1] = Complex(0.0, 0.0);

        for (i = 2; i <= N; i++)
        {
            vectorbeta[N - i] = Cmul(vectoralpha[N - i], Csub(vectorbeta[N - i + 1], Cmul(Complex(0.0, 4.0 / sgorrito), phi[N - i + 1])));
        }

        // Calculo ahora el vector X para tiempo t

        vectorchi[0] = Complex(0.0, 0.0);
        vectorchi[N] = Complex(0.0, 0.0);

        for (i = 1; i < N; i++)
        {
            vectorchi[i] = Cadd(vectorbeta[i - 1], Cmul(vectoralpha[i - 1], vectorchi[i - 1]));
        }

        // Finalmente hallo la nueva onda en tiempo t

        for(i = 0; i <= N; i++)
        {
            phi[i] = Csub(vectorchi[i], phi[i]);
        }

        // Calculamos la probabilidad en el detector a la derecha, y el tiempo para el que se alcanza el primer máximo local

        if (condicion)
        {
            probinteres = 0.0;

            for (i = 4.0 * N / 5.0; i <= N; i++)
            {
                probinteres = probinteres + Cabs(phi[i]) * Cabs(phi[i]);
            }

            if (probinteres < variableauxiliar)
            {
                condicion = false;
                tiempomaximo = t;
                probinteres = variableauxiliar;
            }

            variableauxiliar = probinteres;

            cout << probinteres << endl;
        }

        // Como la probabilidad no cambia ya que no se cambia ninguno de los parámetros iniciales de la onda
        // se ha calculado el valor previo y se ha comparado con un número aleatorio 1000 veces. Así, el coeficiente
        // de transmisión será, tal como cabe esperar, la probabilidad de encontrar la partícula a la derecha del pozo
    }

    // Comparamos la probabilidad obtenida en el intervalo a la derecha con 1000 números aleatorios para obtener
    // el coeficiente de transmisión

    for (i = 1; i <= 1000; i++)
    {
        aleatorio = gsl_rng_uniform(tau);  //número aleatorio real en [0,1]
        
        if (aleatorio < probinteres)
        {
            particulas = particulas + 1;
        }
    }

    transmision = particulas / 1000.0;

    cout << transmision << endl;
    cout << probinteres << endl;
    cout << tiempomaximo << endl;
}