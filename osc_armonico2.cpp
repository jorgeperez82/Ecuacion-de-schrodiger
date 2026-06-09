// Voluntario 3 - Oscilador Armónico Cuántico
#include "./complex.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// funciones iniciales
void CalcularPotencial(vector<double>& V, vector<double>& V_tilde, double w, double dx, int S, double L);
void OndaAutofuncionOG(vector<fcomplex>& phi, int S, double w, double L, double dx);
//funcion mejorada para poder calcular n altos sin reventar
void OndaAutofuncion(vector<fcomplex>& phi, int S, double w, double L, double dx);
void OndaGaussiana(vector<fcomplex>& phi, int S, double dx);

// funcione evolución temporal
void CalcularAlpha(vector<fcomplex>& alpha, vector<double>& V_tilde, int S, double s_tilde);
void CalcularBeta(vector<fcomplex>& beta, vector<fcomplex>& alpha, vector<fcomplex>& phi, vector<double>& V_tilde, double s_tilde, int S);

// funciones observables
double CalcularNorma(vector<fcomplex>& phi, int S, double dx);
double CalcularXMedio(vector<fcomplex>& phi, int S, double dx);
double CalcularXCuadMedio(vector<fcomplex>& phi, int S, double dx);
double IncertidumbreX(double x_med, double x2_med);
double CalcularPMedio(vector<fcomplex>& phi, int S, double dx);
double CalcularPCuadMedio(vector<fcomplex>& phi, int S, double dx);
double IncertidumbreP(double p_med, double p2_med);
double CalcularEMedia(vector<fcomplex>& phi, vector<double>& V, int S, double dx);
double CalcularVMedia(vector<fcomplex>& phi, vector<double>& V, int S, double dx);

int main() {
    double S = 1000.0;
    double dt = 0.0001;
    double w = 200.0;
    double L = 1.0;
    double t = 0.5;
    double n_pasos = t / dt;
    double dx = L/S;

    double s_tilde = dt / (dx * dx);

    vector<double> V(S + 1);
    vector<double> V_tilde(S + 1);
    CalcularPotencial(V, V_tilde, w, dx, S, L);

    vector<fcomplex> phi(S + 1);

    //Descomentar para elegir la funcion inicial 
    ///////////////////////////////////////////
    //OndaAutofuncion(phi, S, w, L, dx);
    OndaGaussiana(phi, S, dx);
    ////////////////////////////////////////

    //Condiciones de contorno
    phi[0] = Complex(0.0, 0.0);
    phi[S] = Complex(0.0, 0.0);

    // Calcular alpha (solo 1 vez)
    vector<fcomplex> alpha(S);
    CalcularAlpha(alpha, V_tilde, S, s_tilde);

    //elegir inteligentemente los nombtes de los ficheros para no compilar 17000 veces
    ofstream fichero1("observables_gaussiana3.txt");
    ofstream fichero2("animacion_gaussiana3.txt");
    //para la energia cinetica y potencil de la gaussiana
    ofstream fichero3("energias_gaussiana3.txt"); 


    //comprobar que se abren los fcheros
    if (!fichero1 || !fichero2 || !fichero3) {
        cerr << "Error al abrir los ficheros." << endl;
        return 1;
    }

    //variables para el calculo de phi con condiciones iniciales
    vector<fcomplex> beta(S);
    vector<fcomplex> chi(S + 1);
    chi[0] = Complex(0.0, 0.0);
    chi[S] = Complex(0.0, 0.0);

    // guaradar instante inicial
    double norma0 = CalcularNorma(phi, S, dx);
    double x_medio0 = CalcularXMedio(phi, S, dx);
    double x2_medio0 = CalcularXCuadMedio(phi, S, dx);
    double p_medio0 = CalcularPMedio(phi, S, dx);
    double p2_medio0 = CalcularPCuadMedio(phi, S, dx);
    double energia0 = CalcularEMedia(phi, V, S, dx);
    double dxdp0 = IncertidumbreX(x_medio0, x2_medio0) * IncertidumbreP(p_medio0, p2_medio0);

    fichero1 << 0.0 << " " << norma0 << " " << x_medio0 << " " << p_medio0 << " " << energia0 << " " << dxdp0 << "\n";

    for (int j = 0; j <= S; j++) {
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        fichero2 << j * dx << " " << phi[j].r << " " << phi[j].i << " " << prob << endl;
    }
    fichero2 << "\n";
    cout << "Instante inicial guardado :)" << endl;

    //guardamos datos de energia potencial y cinetica para la gaussiana
    //podria optimizarse el codigo para no crear un tercer fichero pero ha sido implementado a posteriori ...

    double V_medio0=CalcularVMedia(phi, V, S, dx);
    double T_medio0=CalcularPCuadMedio(phi, S, dx);
    
    fichero3 << 0.0 << " " << T_medio0 << " " << V_medio0 << " " << energia0 << "\n";

    
    
    //evolucion temporal(lo gordo)
    //en muchos casos durante el calculo se opta por hacer ciertas operaciones con complejos sin usar las funciones de la libreria
    //complex.h porque pueden generar algun tipo de imprecision ciertas funciones como Cdiv
    //PD:haciendo debugging se descubrio que la imprecision no era directamente causada por complex.h pero no se volvio a cambiar a
    //funciones de la libreria 
    for (int n = 1; n <= n_pasos; n++) {
        CalcularBeta(beta, alpha, phi, V_tilde, s_tilde, S);

        chi[0] = Complex(0.0, 0.0);

        for (int j = 0; j < S-1; j++) {
            double prod_r = alpha[j].r * chi[j].r - alpha[j].i * chi[j].i;
            double prod_i = alpha[j].r * chi[j].i + alpha[j].i * chi[j].r;
            
            // chi_nueva = prod + beta
            chi[j + 1] = Complex(prod_r + beta[j].r, prod_i + beta[j].i);
        }
        //rrotegemos las condiciones iniciales de algun error de precision computacional volviendo a imponer la condicion inicial
        chi[S]=Complex(0.0,0.0);

        for (int j = 0; j <=S; j++) {
            phi[j] = Complex(chi[j].r - phi[j].r, chi[j].i - phi[j].i);
        }

        phi[0] = Complex(0.0, 0.0);
        phi[S] = Complex(0.0, 0.0);

        // Animación cada 50 pasos
        if (n % 20 == 0 ) {
            for (int j = 0; j <= S; j++) {
                double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
                fichero2 << j * dx << " " << phi[j].r << " " << phi[j].i << " " << prob << endl;
            }
            fichero2 << "\n";
        }

        // Observables (Llamamos a las funciones integradoras)
        double norma = CalcularNorma(phi, S, dx);
        double x_medio = CalcularXMedio(phi, S, dx);
        double x2_medio = CalcularXCuadMedio(phi, S, dx);
        double p_medio = CalcularPMedio(phi, S, dx);
        double p2_medio = CalcularPCuadMedio(phi, S, dx);
        double energia = CalcularEMedia(phi, V, S, dx);
        double dxdp = IncertidumbreX(x_medio, x2_medio) * IncertidumbreP(p_medio, p2_medio);

        fichero1 << n * dt << " " << norma << " " << x_medio << " " << p_medio << " " << energia << " " << dxdp << "\n";
        // para el paquete de ondas gaussiano(editado a posteriori)

        double V_medio=CalcularVMedia(phi, V, S, dx);
        double T_medio=CalcularPCuadMedio(phi, S, dx);
    
        fichero3 << n*dt << " " << T_medio << " " << V_medio << " " << energia << "\n";
    }

    //cerramos ficheros
    fichero1.close();
    fichero2.close();
    fichero3.close();
    cout << "Datos guardados" << endl;
    return 0;
}


void CalcularPotencial(vector<double>& V, vector<double>& V_tilde, double w, double dx, int S, double L) {   
    //funcion que calcula el potencial cuadratico, calculamos V y V_tilde porque se usan los dos durante el codigo
    double x0 = 0.5 * L;
    for (int j = 0; j <= S; j++) {   
        double x = j * dx;
        V[j] = 0.25 * w * w * (x - x0) * (x - x0);
        V_tilde[j] = V[j] * dx * dx;
    }
}

void OndaAutofuncionOG(vector<fcomplex>& phi, int S, double w, double L, double dx) {
    //funcion original, obsoleta porque para n altos se crean errores de precision 
    double norma = 0.0;
    double x0 = 0.5 * L;
    int n = 20; 
    
    for (int j = 0; j <= S; j++) {
        double x = j * dx;
        double expo = -0.25 * w * (x - x0) * (x - x0);
        double arg = sqrt(0.5 * w) * (x - x0);
        double H_n;
        
        if (n == 0) H_n = 1.0;
        else if (n == 1) H_n = 2.0 * arg;
        else if (n == 2) H_n = 4.0 * pow(arg, 2) - 2.0;
        else if (n == 3) H_n = 8.0 * arg * arg * arg - 12.0 * arg;
        //calculamos el polinomio de hermite de orden 20 por recurrencia para que no explote mi pequeño ordenador
        else {
            double h_menos2 = 1.0;       // Equivalente a H_0
            double h_menos1 = 2.0 * arg; // Equivalente a H_1
            H_n = 0.0;
            
            for (int i = 2; i <= n; i++) {
                H_n = 2.0 * arg * h_menos1 - 2.0 * (i - 1) * h_menos2;
                
                // Desplazamos los valores para el siguiente grado
                h_menos2 = h_menos1;
                h_menos1 = H_n;
            }
        }
        
        double phi_real = exp(expo) * H_n;
        phi[j] = Complex(phi_real, 0.0);
        norma += (phi_real * phi_real) * dx;
    }
    
    double cte = 1.0 / sqrt(norma);
    for (int j = 0; j <= S; j++) {
        phi[j] = RCmul(cte, phi[j]);
    }
}

void OndaAutofuncion(vector<fcomplex>& phi, int S, double w, double L, double dx) {
    //funcion para calcular la funcion de onda para ls autofunciones, para la recurrencia en n altos se calcula la
    //funcion de onda por pasos intermedios para evitar errores de precision con potencias tan elevadas
    double norma = 0.0;
    double x0 = 0.5 * L;
    int n = 100; 
    
    for (int j = 0; j <= S; j++) {
        double x = j * dx;
        double expo = -0.25 * w * (x - x0) * (x - x0);
        double arg = sqrt(0.5 * w) * (x - x0);
        double phi_real = 0.0;
        
        if (n == 0) {
            phi_real = exp(expo) * 1.0;
        }
        else if (n == 1) {
            phi_real = exp(expo) * (2.0 * arg);
        }
        else if (n == 2) {
            phi_real = exp(expo) * (4.0 * arg * arg - 2.0);
        }
        else if (n == 3) {
            phi_real = exp(expo) * (8.0 * arg * arg * arg - 12.0 * arg);
        }
        //se amortigua la recurrencia 
        else {
            double phi_menos2 = exp(expo) * 1.0;       // Equivalente a exp * H_0
            double phi_menos1 = exp(expo) * (2.0 * arg); // Equivalente a exp * H_1
            
            for (int i = 2; i <= n; i++) {
                // se multiplica la estructura de Hermite directamente sobre las funciones amortiguadas
                phi_real = 2.0 * arg * phi_menos1 - 2.0 * (i - 1) * phi_menos2;
                
                // filtro de seguridad por si falla todo lo demas 
                if (isnan(phi_real) || isinf(phi_real)) {
                    phi_real = 0.0;
                }
                
                // dsplazamos los valores de las funciones completas
                phi_menos2 = phi_menos1;
                phi_menos1 = phi_real;
            }
        }
        
        // asignamos el valor real limpio de ruido numérico
        phi[j] = Complex(phi_real, 0.0);
        norma += (phi_real * phi_real) * dx;
    }
    
    //normalizamos la funcion de onda SIEMPRE 
    double cte = 1.0 / sqrt(norma);
    for (int j = 0; j <= S; j++) {
        phi[j] = RCmul(cte, phi[j]);
    }
}

void OndaGaussiana(vector<fcomplex>& phi, int S, double dx) {
    //funcion para calcular la funcion de onda gaussiana
    double x0 = 0.5;
    double sigma = 1.0 / 10.0;
    double norma = 0.0;
    
    for (int j = 0; j <= S; j++) {
        double x = j * dx;
        double expo = exp(-pow(x - x0, 2) / (2.0 * sigma * sigma));
        phi[j] = Complex(expo, 0.0);
        norma += (expo * expo) * dx;
    }
    
    //normalizamos
    double cte = 1.0 / sqrt(norma);
    for (int j = 0; j <= S; j++) {
        phi[j] = RCmul(cte, phi[j]);    
    }
}

void CalcularAlpha(vector<fcomplex>& alpha, vector<double>& V_tilde, int S, double s_tilde) {
    //funcion que calcula alpha para poder usarla durante rodo el codigo
    alpha[S - 1] = Complex(0.0, 0.0);
    
    for (int j = S - 1; j > 0; j--) {
        // Denominador = A0 + alpha[j]
        double den_r = -2.0 - V_tilde[j] + alpha[j].r;
        double den_i = 2.0 / s_tilde + alpha[j].i;
        
        // Módulo al cuadrado para dividir a mano
        double mod2 = den_r * den_r + den_i * den_i;
        
        // gamma = 1 / denominador
        double gamma_r = den_r / mod2;
        double gamma_i = -den_i / mod2;
        
        // alpha[j-1] = -gamma 
        alpha[j - 1] = Complex(-gamma_r, -gamma_i);
    }
}

void CalcularBeta(vector<fcomplex>& beta, vector<fcomplex>& alpha, vector<fcomplex>& phi, vector<double>& V_tilde, double s_tilde, int S) {
    //funcion para calcular beta
    beta[S - 1] = Complex(0.0, 0.0);
    
    for (int j = S - 1; j > 0; j--) {
        // Denominador = A0 + alpha[j]
        double den_r = -2.0 - V_tilde[j] + alpha[j].r;
        double den_i = 2.0 / s_tilde + alpha[j].i;
        double mod2 = den_r * den_r + den_i * den_i;
        
        // gamma = 1 / denominador
        double gamma_r = den_r / mod2;
        double gamma_i = -den_i / mod2;
        
        // vector b_j = (4i / s_tilde) * phi[j]
        double b_r = -(4.0 / s_tilde) * phi[j].i;
        double b_i = (4.0 / s_tilde) * phi[j].r;
        
        // suma = b_j - beta[j]
        double suma_r = b_r - beta[j].r;
        double suma_i = b_i - beta[j].i;
        
        // beta[j-1] = suma * gamma (Multiplicación compleja manual)
        double beta_r = suma_r * gamma_r - suma_i * gamma_i;
        double beta_i = suma_r * gamma_i + suma_i * gamma_r;
        
        beta[j - 1] = Complex(beta_r, beta_i);
    }
}

double CalcularNorma(vector<fcomplex>& phi, int S, double dx) {
    //calcula la norma de la funcion(1 para toda funcion normalizada)
    //comprueba que la fisica funciona, 
    double suma = 0.0;
    for (int j = 0; j <= S; j++) {
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        suma += prob * dx;
    }
    return suma;
}

double CalcularXMedio(vector<fcomplex>& phi, int S, double dx) {
    //calculo del valor esperaddo de x
    double suma = 0.0;
    for (int j = 0; j <= S; j++) {
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        suma += (j * dx) * prob * dx;
    }
    return suma;
}

double CalcularXCuadMedio(vector<fcomplex>& phi, int S, double dx) {
    //calculo del valor esperado de x^2, para la incertidumbre de x
    double suma = 0.0;
    for (int j = 0; j <= S; j++) {
        double x = j * dx;
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        suma += (x * x) * prob * dx;
    }
    return suma;
}

double IncertidumbreX(double x_med, double x2_med) {
    //incertidumbre de x
    return sqrt(abs(x2_med - (x_med * x_med)));
}

double CalcularPMedio(vector<fcomplex>& phi, int S, double dx) {
    //calculo del valor esperado del momento
    double suma = 0.0;
    for (int j = 1; j < S; j++) {
        double dR = phi[j + 1].r - phi[j - 1].r;
        double dI = phi[j + 1].i - phi[j - 1].i;
        double valor_p = (phi[j].r * dI - phi[j].i * dR) / (2.0 * dx);
        suma += valor_p * dx;
    }
    return suma;
}

double CalcularPCuadMedio(vector<fcomplex>& phi, int S, double dx) {
    //utilizado para incertidumbre y para el valor esperado de la energia cinetica(iguales por definicion tras el reescalado)
    double suma = 0.0;
    for (int j = 0; j < S-1; j++) {
        double dR = (phi[j + 1].r - phi[j].r) / dx;
        double dI = (phi[j + 1].i - phi[j].i) / dx;
        suma += ((dR * dR) + (dI * dI)) * dx;
    }
    return suma;
}

double IncertidumbreP(double p_med, double p2_med) {
    return sqrt(abs(p2_med - (p_med * p_med)));
}

double CalcularEMedia(vector<fcomplex>& phi, vector<double>& V, int S, double dx) {
    //calculo del valor esperado de la energia E=T+V
    double suma_V = 0.0;
    for (int j = 0; j <= S; j++) {
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        suma_V += V[j] * prob * dx;
    }
    double suma_T = CalcularPCuadMedio(phi, S, dx);
    return suma_T + suma_V;
}

double CalcularVMedia(vector<fcomplex>& phi, vector<double>& V, int S, double dx) {
    double suma_V = 0.0;
    for (int j = 0; j <= S; j++) {
        double prob = (phi[j].r * phi[j].r) + (phi[j].i * phi[j].i);
        suma_V += V[j] * prob * dx;
    }
    return suma_V;
}