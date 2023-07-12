#include "eserc8classes.h"
#include "eserc8functions.h"

using namespace std;

void Reset(int);
void Accumulate(void);
void Averages(int);

int main(int argc, char* argv[]){

    Random rnd;
    initRand(rnd);

    Psi psi(1.,1.);              //wave function and integrand definition

    double inittemp = 1.;
    double temp = inittemp;
    int nstepsTemp = 25;
    int nstepsSA = 100;

    double prob = 0.;          //variables for moves in SA algorithm and Metropolis sampling for integral calculation
    double r = 0.;
    double delta = 7.;           
    double deltaSA = 0.1;
    int acceptedSA = 0; 
    int attemptedSA = 0;
    bool end = false;

    int nblk = 100;
    Metropolis metro(psi, rnd, inittemp, 0., delta, 1E5);
    BlockAverage av(psi, 1, nblk);

    vector<vector<double>> results_old( nblk, vector<double> (4, 0.));      //vector to keep results from the previous iteration

    int wd = 18;                                         //variables for output
    ofstream Annealing, EnergyFinal, Histo, Parameters;

    for(int iblk=1; iblk <= av.nblk; iblk++) {         //calcolo primo integrale

        av.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= metro.nstep; istep++){ 
                    
            metro.SamplingStep();
            av.Measure(metro.xold);
            av.Accumulate();    //Update block averages
                
        }
        av.Averages(iblk);   //Print results for current block
    }
    
    Annealing.open("output_progressive_energy.dat");                    //stampa primo risultato
    Annealing << setw(wd) << av.results[nblk-1][2] << setw(wd) << av.results[nblk-1][3] << endl;

    Parameters.open("output_parameters.dat");   
    Parameters << setw(wd) << temp << setw(wd) << psi.mu << setw(wd) << psi.sigma << endl;

    for(int i = 0; i < nstepsTemp; i++){            //external cicle for temperature

        acceptedSA = 0;
        attemptedSA = 0;

        for(int j = 0; j < nstepsSA; j++){         //cicle for SA algorithm

            results_old = av.results;                       //Salvo i vecchi dati

            av.function.mu += deltaSA * (rnd.Rannyu() - 0.5);
            av.function.sigma += deltaSA * (rnd.Rannyu() - 0.5);     //cambio valori dei parametri in av, per calcolare direttamente l'integrale
            metro.function = av.function;                            //va cambiata anche la psi in metro per il sampling
            
            end = (i == nstepsTemp - 1 and j == nstepsSA - 1);    //bool per non far variare la temperatura all'ultimo step

            for(int iblk=1; iblk <= av.nblk; iblk++) {         //calcolo integrale

                av.Reset(iblk);   //Reset block averages
                for(int istep=1; istep <= metro.nstep; istep++){ 
                    
                    metro.SamplingStep();
                    av.Measure(metro.xold);
                    av.Accumulate();    //Update block averages
                
                }
                av.Averages(iblk);   //Print results for current block
            }

            //ora l'integrale con i nuovi parametri è av.results[nblk-1][2].

            prob = min(1., exp( - (1./temp) * (av.results[nblk-1][2] - results_old[nblk-1][2])));       //calcolo probabilità mossa
            r = rnd.Rannyu();

            if( r < prob ){                     //faccio la mossa se accettata

                psi = av.function;       //parametri vengono aggiornati
                Parameters << setw(wd) << temp << setw(wd) << psi.mu << setw(wd) << psi.sigma << endl;       //stampo parametri se la mossa è stata accettata
                acceptedSA++;

            }else{

                av.function = psi;         //In questo momento psi non è stata cambiata, dunque in psi ci sono i parametri vecchi.
                metro.function = psi; 

                av.results = results_old;       //Elimino l'integrale rifiutato sostituendolo con quello nuovo

            } 
            attemptedSA++;

            Annealing << setw(wd) << av.results[nblk-1][2] << setw(wd) << av.results[nblk-1][3] << endl;      //stampo il risultato della mossa
        }

        std::cout << "--------------------------------------------" << endl;
        std::cout << "Simulation complete at temperature T = " << temp << endl;
        std::cout << "--------------------------------------------" << endl;

        if(!end) {temp *= 0.8;}     //Abbasso la temperatura 
    }

    Annealing.close();
    Parameters.close();

    //Calculation of final integral and histogram of the sampled distribution
    Histo.open("output_histo.dat");

    for(int iblk=1; iblk <= av.nblk; iblk++) {         //calcolo integrale

        av.Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= metro.nstep; istep++){ 
                    
            metro.SamplingStep();
            Histo << metro.xold << endl;
            av.Measure(metro.xold);
            av.Accumulate();    //Update block averages
                
        }
        av.Averages(iblk);   //Print results for current block
    }
    Histo.close();

    EnergyFinal.open("output_energy_final.dat");       //final integral output
    for(int i=0; i < nblk; i++){

        EnergyFinal << setw(wd) << av.results[i][0] <<  setw(wd) << av.results[i][1] << setw(wd) << av.results[i][2] << setw(wd) << av.results[i][3] << endl;
    }
    EnergyFinal.close();

    return 0;
}




