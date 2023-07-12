#include "random.h"
#include "eserc1functions.h"
#include "eserc1classes.h"
#include "Buffonclasses.h"

using namespace std;

int main(int argc, char* argv[]){

    Random rnd;          //inizializzazione numeri casuali
    initRand(rnd);

    //################################# Esercizio 01.1, Punti 1. e 2.

    int nblocks = 300;              //numero di blocchi
    int nthrows_1 = 1000;              //numero di lanci per blocco

    r fun;

    BlockAverage av(fun, 2, nblocks);

    for(int iblk=1; iblk <= av.nblk; iblk++) {
    
        av.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nthrows_1; istep++){ 
            
            av.Measure(rnd.Rannyu());
            av.Accumulate();    //Update block averages
          
        }
        av.Averages(iblk);   //Print results for current block
    }

    //############################################

    //################################## Esercizio 01.2, Chi quadro

    int nbins = 100;              //numero di cestini
    int nthrows = 1E4;            //numero di lanci
    int num = 1E4;                //numero di chi2 calcolati

    ofstream chiout;
    chiout.open("chiout.dat");

    double r = 0.;
    double chi2_ = 0.;
    vector<int> conteggi(nbins);

    for(int i=0; i<num; i++){
            
        chi2_ = 0.;
        fill(conteggi.begin(),conteggi.end(), 0);  //vettore per i conteggi dei numeri in ciascun bin

        for(int j=0; j<nthrows; j++){     

            r = rnd.Rannyu();  
            conteggi[ floor(r * (double) nbins) ]++;      //aumenta contatore del bin
        }
        for(int j=0; j<nbins; j++){

            chi2_ += pow((conteggi[j] - nthrows/nbins),2)/(nthrows/nbins);        //calcolo chi2
        }

        chiout << setw(15) << chi2_ << endl;
    }
    chiout.close();

    //#################################################################

    //######################################### Esercizio 01.2

    int nthrows_2 = 10000;              //numero di lanci 

    vector<double> countsStandard {0,0,0,0};      
    vector<double> countsExp {0,0,0,0};
    vector<double> countsLorentz {0,0,0,0};

    ofstream Histo_standard("Histo_standard.dat");
    ofstream Histo_exp("Histo_exp.dat");
    ofstream Histo_lorentz("Histo_lorentz.dat");

    int wd = 18;

    for(int i = 0; i < nthrows_2; i++){

        countsStandard[0] += rnd.Rannyu(0.,1.);         //N = 1
        countsExp[0] += rnd.Exp(0.,1.);
        countsLorentz[0] += rnd.Lorentz(0.,1.);

        for(int j = 0; j < 2; j++){                   // N = 2

            countsStandard[1] += rnd.Rannyu(0.,1.);
            countsExp[1] += rnd.Exp(0.,1.);
            countsLorentz[1] += rnd.Lorentz(0.,1.); 
        }

        for(int j = 0; j < 10; j++){                   // N = 10

            countsStandard[2] += rnd.Rannyu(0.,1.);
            countsExp[2] += rnd.Exp(0.,1.);
            countsLorentz[2] += rnd.Lorentz(0.,1.); 
        }

        for(int j = 0; j < 100; j++){                   // N = 100

            countsStandard[3] += rnd.Rannyu(0.,1.);
            countsExp[3] += rnd.Exp(0.,1.);
            countsLorentz[3] += rnd.Lorentz(0.,1.); 
        }

        Histo_standard << setw(wd) << countsStandard[0] << setw(wd) << countsStandard[1] / 2. << setw(wd) << countsStandard[2] / 10. << setw(wd) << countsStandard[3] / 100. << endl;
        Histo_exp << setw(wd) << countsExp[0] << setw(wd) << countsExp[1] / 2. << setw(wd) << countsExp[2] / 10. << setw(wd) << countsExp[3] / 100. << endl;
        Histo_lorentz << setw(wd) << countsLorentz[0] << setw(wd) << countsLorentz[1] / 2. << setw(wd) << countsLorentz[2] / 10. << setw(wd) << countsLorentz[3] / 100. << endl;

        fill(countsStandard.begin(), countsStandard.end(), 0.);
        fill(countsExp.begin(), countsExp.end(), 0.);
        fill(countsLorentz.begin(), countsLorentz.end(), 0.);

    }

    Histo_standard.close();
    Histo_exp.close();
    Histo_lorentz.close();

    //#########################################################

    //################################################ Esercizio 01.3, Buffon

    int nblockspi = 10000;              //numero di blocchi
    int nthrowspi = 1E5;              //numero di lanci per blocco

    pi buffon(rnd, 0.75, 1.);

    BlockAveragePi ave(buffon, nblockspi, 1);

    for(int iblk=1; iblk <= nblockspi; iblk++) {
    
        ave.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nthrowspi; istep++){ 

            ave.Measure();
            ave.Accumulate();    //Update block averages
          
        }
        ave.Averages(iblk);   //Print results for current block
    }

    return 0;
}