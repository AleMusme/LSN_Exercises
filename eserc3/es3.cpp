#include "random.h"
#include "eserc3functions.h"
#include "eserc3classes.h"

using namespace std;

int main(int argc, char* argv[]){

    Random rnd;          //inizializzazione numeri casuali
    initRand(rnd);

    Options opt(rnd);

    int nblocks = 100;              //numero di blocchi
    int nthrows = opt.GetNum() / nblocks;   //numero di lanci per blocco

    BlockAverage av(opt, 4, nblocks); 

    for(int iblk=1; iblk <= av.nblk; iblk++) {
    
        av.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nthrows; istep++){ 
            
            av.Measure( (iblk - 1) * nthrows + istep);
            av.Accumulate();    //Update block averages
          
        }
        av.Averages(iblk);   //Print results for current block
    }

    return 0;
}