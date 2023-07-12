#ifndef __eserc10functions__
#define __eserc10functions__

#include "mpi.h"
#include "random.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void initRand(Random &rnd, int plus);
double error(double av1, double av2, int n);
int Pbc(vector<int> v, int r);
int Pbcnew(vector<int> &v, int r);

void initRand(Random &rnd, int plus){

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
    seed[3] += plus;

    rnd.SetRandom(seed,p1,p2);
    Seed.close();

}

void Shuffle(Random& rnd, vector<int> &v){

    int k = 0;
    int l = 0;

    while( k==l ){
        k = floor(rnd.Rannyu(1., v.size()));
        l = floor(rnd.Rannyu(1., v.size()));
    }
    swap(v[k], v[l]);
}

int Pbcnew(vector<int> &v, int r){

    return r - v.size() * floor( (double) r/v.size());
}

double error(double sum, double sum2, int iblk){  //funzione per il calcolo delle incertezze 

    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

int Pbc(vector<int> v, int r){   //r is an index. the function gives an index, which is r, cycling back to element 1 if index exceeds length.

    return r - (v.size() - 1) * floor( (double) r/v.size() );
}

#endif