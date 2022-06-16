#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <fstream>
#include <cstring>

using namespace std;

ofstream pisz;

//Number of dimensions
#define IL_ZMIENNYCH 30
//Number of generations
#define MAX_GENERATION 10000
//PI number
#define M_PI 3.14159265358979323846
//Population size
#define POPSIZE 20
//Number of descendants
#define LPOTOMKOW 140
//History size
#define HISTORIA 5
//Number of rules
#define MAX_REGUL 4

//Maximum fitness value
#define MAX 1000

//Generation number
int generation;
//mean value for mutation
auto mean=0;
//stddev for mutation
auto stddev = 0.5;
double stddev_tmp;

int liczba_mutacji_pozytywnych[HISTORIA];
int il_mut;
double historia_fitness[HISTORIA];
double avg_fitness;


class genotype
{
public:
  double gene[IL_ZMIENNYCH];		//gene string
  double fitness;					//individual's fitness
  double oblicz_dopasowanie();      //method for calculating the fitness
};

genotype osobnik[POPSIZE];
genotype potomek[LPOTOMKOW+POPSIZE];

/**************************************************************************/
// memberships for fuzzy logic
/**************************************************************************/
double lmp5_procent;
double avg_fitness_ratio;

double mf_liczba_mutacji[2];
double mf_wsp_przystosowania[2];
double mf_wynik[2];

double f_mut_m, f_mut_d;
double f_przyst_m, f_przyst_d;

double r[MAX_REGUL];

/**************************************************************************/
// method for calculating the fitness
/**************************************************************************/

double genotype::oblicz_dopasowanie()
{
    int k;
	double dopasowanie = 10 * IL_ZMIENNYCH;
    for (k=0; k<IL_ZMIENNYCH; k++)
    {
        dopasowanie = dopasowanie + gene[k]*gene[k]-10*cos(2*M_PI*gene[k]);
    }

return dopasowanie;
};

/**************************************************************************/
//Required C++11
//Settings/Compiler/Have g++ follow C++11
/**************************************************************************/
double losuj(double m, double s)
{
    double wynik;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator(seed);
    normal_distribution<double> distribution(m, s);

    wynik = distribution(generator);
    return wynik;
}

/**************************************************************************/
void initialize(void)
{
 int i, mem;

    for (mem=0; mem<POPSIZE; mem++)
    	{
        for (i=0; i<IL_ZMIENNYCH; i++)
            {
                osobnik[mem].gene[i]=(double)(rand()%1000)/1000*10.24 - 5.12;
            }
        osobnik[mem].fitness = MAX - osobnik[mem].oblicz_dopasowanie();
    	}


    for (mem=0; mem<POPSIZE; mem++)
    	{
        cout << "individual " << mem << endl;
        for (i=0; i<IL_ZMIENNYCH; i++)
            {
                cout << osobnik[mem].gene[i] << ", ";
            }
        cout << endl;
        cout << osobnik[mem].fitness << endl;
    	}


    avg_fitness=0;
    for (i=0; i<POPSIZE; i++)
    {
        avg_fitness += osobnik[i].fitness;
    }
    avg_fitness = avg_fitness/POPSIZE;
    historia_fitness[0] = avg_fitness;
};

/**************************************************************************/
void generuj_nowego()
{
    int i, mem, rodzic;
    int gen_rand;

    il_mut=0;

    for (mem=0; mem<LPOTOMKOW; mem++)
    	{
        gen_rand = rand()%IL_ZMIENNYCH;
        rodzic = rand()%POPSIZE;

        for (i=0; i<IL_ZMIENNYCH; i++)
            {
                potomek[mem].gene[i] = osobnik[rodzic].gene[i];
            }
        potomek[mem].gene[gen_rand] = potomek[mem].gene[gen_rand] + losuj(mean, stddev_tmp);
        if (potomek[mem].gene[gen_rand]<-5.12) {potomek[mem].gene[gen_rand]=-5.12;};
        if (potomek[mem].gene[gen_rand]>5.12) {potomek[mem].gene[gen_rand]=5.12;};

        potomek[mem].fitness = MAX - potomek[mem].oblicz_dopasowanie();

        if (potomek[mem].fitness > osobnik[rodzic].fitness)
            {
                il_mut++;
            }
    	}

    for (mem=0; mem<POPSIZE; mem++)
    	{
        for (i=0; i<IL_ZMIENNYCH; i++)
            {
                potomek[LPOTOMKOW+mem].gene[i] = osobnik[mem].gene[i];
            }
        potomek[LPOTOMKOW+mem].fitness = osobnik[mem].fitness;
    	}
}

/*************************************************************************/
// Selection a new population
/*************************************************************************/

void select()
{

    int i, mem;
    int poz_max = 0;
    double max = potomek[0].fitness;
    avg_fitness=0;

    for (mem=0; mem<POPSIZE; mem++)
    {
        poz_max = 0;
        max = potomek[0].fitness;
        for (i=1; i<LPOTOMKOW+POPSIZE; i++)
        {
           if (potomek[i].fitness >= max)
           {
            poz_max = i;
            max = potomek[i].fitness;
           }
        }
        for (i=0; i<IL_ZMIENNYCH; i++)
        {
           osobnik[mem].gene[i] = potomek[poz_max].gene[i];
        }
        osobnik[mem].fitness = potomek[poz_max].fitness;
        potomek[poz_max].fitness = 0;

        avg_fitness += osobnik[mem].fitness;
    }
    avg_fitness = avg_fitness/POPSIZE;

}
/**************************************************************************/
void mf()
{
    mf_liczba_mutacji[0] = 0.01;
    mf_liczba_mutacji[1] = 5;

    mf_wsp_przystosowania[0] = 1.0;
    mf_wsp_przystosowania[1] = 1.1;

    mf_wynik[0] = 0.95;
    mf_wynik[1] = 1.05;

}
/**************************************************************************/
void rozmyj(double mut, double przyst)
{
    if (mut<mf_liczba_mutacji[0])
    {
        f_mut_m=1;
        f_mut_d=0;
    }
    else
    {
        if (mut<mf_liczba_mutacji[1])
        {
            f_mut_m=(mut-mf_liczba_mutacji[0])/(mf_liczba_mutacji[1]-mf_liczba_mutacji[0]);
            f_mut_d=1-f_mut_m;
        }
        else
        {
            f_mut_m=0;
            f_mut_d=1;
        }
    }

    if (przyst<mf_wsp_przystosowania[0])
    {
        f_przyst_m=1;
        f_przyst_d=0;
    }
    else
    {
        if (przyst<mf_wsp_przystosowania[1])
        {
            f_przyst_m=(przyst-mf_wsp_przystosowania[0])/(mf_wsp_przystosowania[1]-mf_wsp_przystosowania[0]);
            f_przyst_d=1-f_przyst_m;
        }
        else
        {
            f_przyst_m=0;
            f_przyst_d=1;
        }
    }

}
/**************************************************************************/
void wnioskowanie()
{
    r[0] = f_mut_m;
    r[1] = f_mut_d;
    r[2] = f_przyst_m;
    r[3] = f_przyst_d;
}

/**************************************************************************/
double wyostrz()
{
    double m, d, wynik;
    m = r[0] + r[2];
    d = r[1] + r[3];
    wynik = d/(m+d);
    wynik = mf_wynik[0] + wynik * (mf_wynik[1]-mf_wynik[0]);

    return wynik;
}

/**************************************************************************/
int main()
{
    int i;
    int lmp5;
    double lmp5_procent_min = 10000, lmp5_procent_max= 0;
    double avg_fitness5;
    double avg_fitness_ratio_min = 10000, avg_fitness_ratio_max = 0;

    srand(time(0));

    pisz.open("log.txt");
      if(!pisz)
        {cout << "Input file error" << endl;
        return 1;
        }

    stddev_tmp = stddev;
    mf();

    generation = 0;
    initialize();

    do
    {
    generuj_nowego();
    select();

/**************************************************************************/
//Fitness history analysis
/**************************************************************************/

    for (i=HISTORIA-2; i>=0; i--)
    {
        historia_fitness[i+1] = historia_fitness[i];
    }
    historia_fitness[0] = avg_fitness;

    if (generation > 4)
    {
        avg_fitness5=0;
        for (i=0; i<HISTORIA; i++)
            {
                avg_fitness5 += historia_fitness[i];
            }
        avg_fitness5 = avg_fitness5/HISTORIA;

        avg_fitness_ratio = historia_fitness[0] / avg_fitness5;

        if (avg_fitness_ratio < avg_fitness_ratio_min)
        {
            avg_fitness_ratio_min = avg_fitness_ratio;
        }

        if (avg_fitness_ratio > avg_fitness_ratio_max)
        {
            avg_fitness_ratio_max = avg_fitness_ratio;
        }
    }

/**************************************************************************/
//Analysis of the number of positive mutations
/**************************************************************************/

    for (i=HISTORIA-2; i>=0; i--)
    {
        liczba_mutacji_pozytywnych[i+1] = liczba_mutacji_pozytywnych[i];
    }
    liczba_mutacji_pozytywnych[0] = il_mut;


    if (generation > 4)
    {
        lmp5=0;
        for (i=0; i<HISTORIA; i++)
        {
            lmp5 += liczba_mutacji_pozytywnych[i];
        }
        lmp5_procent = 100*(double)(lmp5)/HISTORIA/LPOTOMKOW;

        if (lmp5_procent < lmp5_procent_min)
        {
            lmp5_procent_min = lmp5_procent;
        }

        if (lmp5_procent > lmp5_procent_max)
        {
            lmp5_procent_max = lmp5_procent;
        }
    }

/**************************************************************************/
//Fuzzy Logic analysis
/**************************************************************************/
    rozmyj(lmp5_procent, avg_fitness_ratio);
    wnioskowanie();

    stddev_tmp = stddev * wyostrz();
    cout << "stddev " << stddev << endl;

    generation++;
    cout << "generation " << generation << " fitness " << MAX - osobnik[0].fitness << endl;
    pisz << generation << ". " << MAX-osobnik[0].fitness << endl;
    }while(generation < MAX_GENERATION);
    ///}while(osobnik[0].fitness <= MAX - 0.001);

        cout << "The best - genes" << endl;
        for (i=0; i<IL_ZMIENNYCH; i++)
        {
            cout << osobnik[0].gene[i] << ", ";
        }
        cout << endl;
        cout << "fitness " << MAX - osobnik[0].fitness << endl;


pisz.close();
    return 0;
}






