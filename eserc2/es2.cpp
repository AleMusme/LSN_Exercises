#include "eserc2functions.h"
#include "eserc2classes.h"
#include "RWclasses.h"
#include "random.h"

using namespace std;

void initRand(Random &rnd);
double error(double av1, double av2, int n);

int main(int argc, char* argv[]){

    //################################# Esercizio 02.1

    int nblocks = 1000;              //numero di blocchi
    int nthrows = 1000;              //numero di lanci per blocco

    Random rnd;          //inizializzazione numeri casuali
    initRand(rnd);

    f fun;

    BlockAverage av(fun, 2, nblocks);

    for(int iblk=1; iblk <= av.nblk; iblk++) {
    
        av.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nthrows; istep++){ 

            av.Measure(rnd.Rannyu(), rnd.distr());
            av.Accumulate();    //Update block averages
        }
        av.Averages(iblk);   //Print results for current block
    }    

    //##################################################

    //################################# Esercizio 02.2

    int nwalkers = 10000;  //numero di random walk totali
    int nsteps = 100;  //numero di step del random walk
    nblocks = 100;  //numero di blocchi
    int nvalues = nwalkers/nblocks;   //numero di random walkers in ogni blocco

    path dis(rnd, nwalkers);
    ofstream Discreto("output_discreto.dat");
    BlockAverageRW avDis(dis, 2, nblocks);

    path cont(rnd, nwalkers);
    ofstream Continuo("output_continuo.dat");
    BlockAverageRW avCont(cont, 2, nblocks);

    for(int imove=0; imove < nsteps; imove++){

        for(int iblk=1; iblk <= avDis.nblk; iblk++) {
    
            avDis.Reset(iblk);   //Reset block averages
            avCont.Reset(iblk);
            for(int istep=1; istep <= nvalues; istep++){ 

                avDis.Measure((iblk - 1)* nvalues + istep);
                avDis.Accumulate();    //Update block averages

                avCont.Measure((iblk - 1)* nvalues + istep);
                avCont.Accumulate();    
            }
            avDis.Averages(iblk, Discreto);   //Print results for final block
            avCont.Averages(iblk, Continuo);
        }

        dis.StepDiscreto();       //tutti i walkers avanzano di uno step
        avDis.p = dis;

        cont.StepContinuo();
        avCont.p = cont;
    }

    Discreto.close();
    Continuo.close();

    //##################################################

    return 0;
}