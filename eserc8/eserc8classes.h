#ifndef __eserc8classes__
#define __eserc8classes__

#include "random.h"
#include "eserc8functions.h"

using namespace std;

class Psi {

    public:

        Psi(double s = 1., double em = 1.) {
            sigma = s;
            mu = em;
            h = 1.;
            m = 1.;
        }
        
        double Eval(double x){
            return exp( - pow((x - mu), 2)/(2*sigma*sigma) ) + exp( - pow((x + mu), 2 )/(2*sigma*sigma) );
        }

        double Evalmod2(double x){
            return pow( (exp( - pow((x - mu), 2)/(2*sigma*sigma) ) + exp( - pow((x + mu), 2 )/(2*sigma*sigma) ) )  , 2);
        }

        double EvalIntegrand(double x){
            return h*h/(2*m*pow(sigma,4)) * (sigma*sigma - mu*mu - x*x + 2*x*mu*tanh(x*mu/(sigma*sigma))) +
                (pow(x,4) - 5./2.*x*x);  
        }

    public:
        double sigma, mu;
        double h; 
        double m;
        double beta;

};


class Metropolis{

    public:

        Metropolis(){}

        Metropolis(Psi& p, Random& r, double t = 10., double x0 = 0., double d = 1., int n = 1E5) :
            nstep{n}, xn{x0}, xold{x0}, delta{d}, temp{t}, function{p}, rnd{r} {

            accepted = 0;
            attempted = 0;
            acceptedSA = 0;
            attemptedSA = 0;
        }

        ~Metropolis(){}

        void SamplingStep(){

            xn = xold + delta*(rnd.Rannyu() - 0.5);

            double A = min(1. , function.Evalmod2(xn)/function.Evalmod2(xold));

            double r = rnd.Rannyu();

            if( r < A ){

                xold = xn;
                accepted++;
            }
            attempted++;
        }

    public:
        int nstep;
        int accepted, attempted;
        double xn, xold;
        double delta;
        double temp;
        double deltaSA;

        int acceptedSA, attemptedSA;

        Psi function;
        Random rnd;

};

class BlockAverage{

    public:

        BlockAverage(Psi& p, int n = 1, int b = 20) :
            results(b, vector<double> (4, 0.)), n_props{n}, nblk{b}, blk_av(n_props), glob_av(n_props), glob_av2(n_props), walker(n_props), function{p} 
            {
                blk_norm = 0.;
                stima = 0.;
                errore = 0.;
                result = 0.;
                result_error = 0.;
                end = false;
            }

        void Reset(int iblk){ //Reset block averages
    
            if(iblk == 1){

                for(int i=0; i<n_props; ++i)
                {
                    glob_av[i] = 0;
                    glob_av2[i] = 0;
                }
            }

            for(int i=0; i<n_props; ++i){

                blk_av[i] = 0;
            }

            blk_norm = 0;
            stima = 0.;
            errore = 0.;
            result = 0.;
            result_error = 0.;
        }

        void Measure(double x){

            walker[ie] = function.EvalIntegrand(x);
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block

            stima = blk_av[ie]/blk_norm; //Energy
            glob_av[ie] += stima;
            glob_av2[ie] += stima*stima;
            errore = error(glob_av[ie],glob_av2[ie],iblk);

            results[iblk-1][0] = iblk;
            results[iblk-1][1] = stima;
            results[iblk-1][2] = glob_av[ie] / iblk;
            results[iblk-1][3] = errore;

        }

    public:
    
        int n_props, nblk;
        const int ie=0;

        vector<double>  blk_av, glob_av, glob_av2;
        vector<double> walker;
        vector <vector<double>> results;
        double stima;
        double errore;
        double blk_norm;
        double result;
        double result_error;
        bool end;

        Psi function;

};


#endif 