#ifndef __eserc3functions__
#define __eserc3functions__

#include "random.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void initRand(Random &rnd);
double error(double av1, double av2, int n);

void initRand(Random &rnd){

    ifstream Primes("Primes");

    //Read seed for random numbers
    int p1, p2;
    Primes >> p1 >> p2 ;
    Primes.close();

    //if(restart) Seed.open("seed.out");
    //else Seed.open("seed.in");

    int seed[4];

    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

}

double error(double sum, double sum2, int iblk){  //funzione per il calcolo delle incertezze 

    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}






#endif