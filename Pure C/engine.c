/* **********************************************************************
   This is the low-level implementation of Schrodinger equation solver
   using C with Runge-Kutta forth-order method. The main purpose of
   this program is to find out the highest possible efficiency to solve
   Schrodinger equation with Runge-Kutta method.

       Author: Zhang Chang-kai
       E-mail: phy.zhangck@gmail.com

   This program is license under GPL-3.0 License.
********************************************************************** */

#include "stdio.h"
#include "math.h"

/* basic parameters */
#define TIME_STEP 4e-6
#define SPACE_STEP 0.02
#define TIME_STEPS 100000
#define SPACE_STEPS 51
#define TARGET 9997

/* struct for complex numbers */
struct complex
{
    long double real;
    long double imag;
};

/* calculate the central difference */
long double central(long double l, long double c, long double r)
{
    long double df;
    df = l + r - 2 * c;
    df = df / pow(SPACE_STEP, 2);
    return df;
}

/* potential */
long double potential(void)
{
    long double V = 0.0;
    return V;
}

/* calculate the hamiltonian */
long double hamiltonian(long double l, long double c, long double r)
{
    long double ham;
    ham = (-1) * central(l, c, r) + potential() * c;
    return ham;
}

/* Runge-Kutta forth-order estimation */
struct complex* estimate(struct complex psi0[SPACE_STEPS])
{
    /* k1, k2, k3, k4, k, psi */
    struct complex k1[SPACE_STEPS];
    struct complex k2[SPACE_STEPS];
    struct complex k3[SPACE_STEPS];
    struct complex k4[SPACE_STEPS];
    struct complex k[SPACE_STEPS];
    struct complex psi[SPACE_STEPS];

    /* time step size and number of space steps */
    long double h = TIME_STEP;
    int nx = SPACE_STEPS;

    /* calculate k1 */
    for (int i = 1; i < nx - 1; i++)
    {
        k1[i].real = hamiltonian(psi0[i-1].imag,
                                 psi0[i].imag,
                                 psi0[i+1].imag);
        k1[i].imag = hamiltonian(psi0[i-1].real,
                                 psi0[i].real,
                                 psi0[i+1].real) * (-1);
    }
    k1[0].real = k1[0].imag = k1[nx-1].real = k1[nx-1].imag = 0.0;

    /* calculate k2 */
    for (int i = 1; i < nx - 1; i++)
    {
        k2[i].real = hamiltonian(psi0[i-1].imag + k1[i-1].imag*h/2,
                                 psi0[i].imag + k1[i].imag*h/2,
                                 psi0[i+1].imag + k1[i+1].imag*h/2);
        k2[i].imag = hamiltonian(psi0[i-1].real + k1[i-1].real*h/2,
                                 psi0[i].real + k1[i].real*h/2,
                                 psi0[i+1].real + k1[i+1].real*h/2) * (-1);
    }
    k2[0].real = k2[0].imag = k2[nx-1].real = k2[nx-1].imag = 0.0;

    /* calculate k3 */
    for (int i = 1; i < nx - 1; i++)
    {
        k3[i].real = hamiltonian(psi0[i-1].imag + k2[i-1].imag*h/2,
                                 psi0[i].imag + k2[i].imag*h/2,
                                 psi0[i+1].imag + k2[i+1].imag*h/2);
        k3[i].imag = hamiltonian(psi0[i-1].real + k2[i-1].real*h/2,
                                 psi0[i].real + k2[i].real*h/2,
                                 psi0[i+1].real + k2[i+1].real*h/2) * (-1);
    }
    k3[0].real = k3[0].imag = k3[nx-1].real = k3[nx-1].imag = 0.0;

    /* calculate k4 */
    for (int i = 1; i < nx - 1; i++)
    {
        k4[i].real = hamiltonian(psi0[i-1].imag + k3[i-1].imag*h/2,
                                 psi0[i].imag + k3[i].imag*h/2,
                                 psi0[i+1].imag + k3[i+1].imag*h/2);
        k4[i].imag = hamiltonian(psi0[i-1].real + k3[i-1].real*h/2,
                                 psi0[i].real + k3[i].real*h/2,
                                 psi0[i+1].real + k3[i+1].real*h/2) * (-1);
    }
    k4[0].real = k4[0].imag = k4[nx-1].real = k4[nx-1].imag = 0.0;

    /* calculate k = (k1 + 2*k2 + 2*k3 + k4) / 6 */
    for (int i = 0; i < nx; i++)
    {
        k[i].real = k1[i].real + 2*k2[i].real + 2*k3[i].real + k4[i].real;
        k[i].imag = k1[i].imag + 2*k2[i].imag + 2*k3[i].imag + k4[i].imag;
        k[i].real /= 6; k[i].imag /= 6;
    }

    /* estimate psi */
    for (int i = 0; i < nx; i++)
    {
        psi[i].real = psi0[i].real + k[i].real*h;
        psi[i].imag = psi0[i].imag + k[i].imag*h;
    }

    /* return array pointer */
    return psi;
}

int main(void)
{
    /* store simulation result */
    static struct complex psi[10000][SPACE_STEPS];
    static long double psir[TIME_STEPS][SPACE_STEPS];
    static long double psii[TIME_STEPS][SPACE_STEPS];

    /* initial condition */
    for (int i = 0; i < SPACE_STEPS; i++)
    {
        psi[0][i].real = exp(-pow(i*SPACE_STEP - 0.5, 2) / 0.01);
        psi[0][i].imag = 0.0;
    }

    /* update data */
    /* the "mod 10000" is necessary for macOS since there
       is a limit for stack size */
    for (int i = 0; i < TIME_STEPS; i++)
    {
        struct complex* reslt = estimate(psi[i%10000]);
        for (int j = 0; j < SPACE_STEPS; j++)
        {
            psi[(i+1)%10000][j] = *(reslt+j);
            psir[i][j] = (*(reslt+j)).real;
            psii[i][j] = (*(reslt+j)).imag;
        }
    }

    /* output real part */
    for (int j = 0; j < TIME_STEPS; j+=100) {
        printf("%s", "[");
        for (int i = 0; i < SPACE_STEPS - 1; i++)
            { printf("%.18Lf, ", psir[j][i]); }
        printf("%.18Lf", psir[j][SPACE_STEPS-1]);
        printf("%s", "]\n");
    }

    /* print mark */
    printf("%s\n", "switch");

    /* output imaginary part */
    for (int j = 0; j < TIME_STEPS; j+=100) {
        printf("%s", "[");
        for (int i = 0; i < SPACE_STEPS - 1; i++)
            { printf("%.18Lf, ", psii[j][i]); }
        printf("%.18Lf", psii[j][SPACE_STEPS-1]);
        printf("%s", "]\n");
    }

    return 0;
}
