#ifndef __eserc3classes__
#define __eserc3classes__

#include "random.h"
#include "eserc3functions.h"

using namespace std;

class Options {

    public:

        Options(Random& rand, int n=1E4, double T_v=1., double S0_v=100., double r_v=0.1, double K_v=100., double sigma_v=0.25) :
            callDirect(n,0), putDirect(n,0), callDiscrete(n,0), putDiscrete(n,0), num{n}, T{T_v}, S0{S0_v}, r{r_v}, K{K_v}, sigma{sigma_v}, rnd{rand} {

                double S_T = 0.;  //sampling of S(T)

                double sumCallDirect = 0.;
                double sumPutDirect = 0.;

                int nsteps = 100;           //parameters for discrete sampling
                double interval = T / (double) nsteps;
                double S_dis = S0;
                double sumCallDiscrete = 0.;
                double sumPutDiscrete = 0.;

                double nthrows = 1000;       //number of throws to calculate S_0

                for(int i=0; i < num; i++){          

                    S_T = 0.;
                    sumCallDirect = 0.;
                    sumPutDirect = 0.;

                    sumCallDiscrete = 0.;
                    sumPutDiscrete = 0.;

                    for(int j=0; j < nthrows; j++){

                        S_T = S0 * exp( (r - 0.5*sigma*sigma)* T + sigma*rnd.Gauss(0., T)*sqrt(T) );   //Price at time T

                        S_dis = S0;

                        for(int k=0; k < nsteps; k++){
                            
                            S_dis = S_dis * exp( (r - 0.5*sigma*sigma) * interval + sigma * rnd.Gauss(0.,1.) * sqrt(interval) );       //recursive formula for price at discrete times
                        }

                        sumCallDirect += exp(-r*T) * max(0., S_T - K);     //call option, Direct Sampling
                        sumPutDirect += exp(-r*T) * max(0., K - S_T);     //put option, Direct Sampling

                        sumCallDiscrete += exp(-r*T) * max(0., S_dis - K);  //call option, Discrete Sampling
                        sumPutDiscrete += exp(-r*T) * max(0., K - S_dis);     //put option, Discrete Sampling

                    }

                    callDirect[i] = sumCallDirect / (double) nthrows;
                    putDirect[i] = sumPutDirect / (double) nthrows;

                    callDiscrete[i] = sumCallDiscrete / (double) nthrows;
                    putDiscrete[i] = sumPutDiscrete / (double) nthrows;
                }
            }

        int GetNum(){
            return num;
        }
        vector<double> GetCallDir(){
            return callDirect;
        }
        vector<double> GetPutDir(){
            return putDirect;
        }
        vector<double> GetCallDis(){
            return callDiscrete;
        }
        vector<double> GetPutDis(){
            return putDiscrete;
        }

        

    private:                             
        vector<double> callDirect;
        vector<double> putDirect;

        vector<double> callDiscrete;
        vector<double> putDiscrete;

        int num;

        double T;
        double S0;
        double r;
        double K;
        double sigma;

        Random rnd;

};

class BlockAverage{

    public:

        BlockAverage(Options& fun, int n = 4, int b = 20):
            n_props{n}, nblk{b}, blk_av(n_props), glob_av(n_props), blk_av2(n_props), glob_av2(n_props), walker(n_props), function{fun} {
                blk_norm = 0.;
                stima = 0.;
                errore = 0.;
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
        }

        void Measure(int i){

            walker[icdr] = function.GetCallDir()[i];
            walker[ipdr] = function.GetPutDir()[i];
            walker[icds] = function.GetCallDis()[i];
            walker[ipds] = function.GetPutDis()[i];
        }

        void Accumulate(void){ //Update block averages
 
            for(int i=0; i<n_props; ++i)
            {
                blk_av[i] = blk_av[i] + walker[i];
            }
            blk_norm = blk_norm + 1.0;
            
        }

        void Averages(int iblk){ //Print results for current block
        
            ofstream CallDirect, PutDirect, CallDiscrete, PutDiscrete; 
            const int wd=12;
            
            //cout << "Block number " << iblk << endl;

            CallDirect.open("output_CallDirect.dat",ios::app);
            PutDirect.open("output_PutDirect.dat",ios::app);
            CallDiscrete.open("output_CallDiscrete.dat",ios::app);
            PutDiscrete.open("output_PutDiscrete.dat",ios::app);

            stima = blk_av[icdr]/blk_norm; // CallDirect
            glob_av[icdr] += stima;
            glob_av2[icdr] += stima*stima;
            errore = error(glob_av[icdr],glob_av2[icdr],iblk);

            CallDirect << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[icdr]/(double)iblk << setw(wd) << errore << endl;

            stima = blk_av[ipdr]/blk_norm; // PutDirect
            glob_av[ipdr] += stima;
            glob_av2[ipdr] += stima*stima;
            errore = error(glob_av[ipdr],glob_av2[ipdr],iblk);

            PutDirect << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[ipdr]/(double)iblk << setw(wd) << errore << endl;

            stima = blk_av[icds]/blk_norm; // CallDiscrete
            glob_av[icds] += stima;
            glob_av2[icds] += stima*stima;
            errore = error(glob_av[icds],glob_av2[icds],iblk);

            CallDiscrete << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[icds]/(double)iblk << setw(wd) << errore << endl;

            stima = blk_av[ipds]/blk_norm; // PutDiscrete
            glob_av[ipds] += stima;
            glob_av2[ipds] += stima*stima;
            errore = error(glob_av[ipds],glob_av2[ipds],iblk);

            PutDiscrete << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[ipds]/(double)iblk << setw(wd) << errore << endl;

            //cout << "----------------------------" << endl << endl;

            CallDirect.close();
            PutDirect.close();
            CallDiscrete.close();
            PutDiscrete.close();

        }


    public:
    
        int n_props, nblk;
        const int icdr=0;
        const int ipdr=1;
        const int icds=2;
        const int ipds=3;

        vector<double>  blk_av, glob_av, blk_av2, glob_av2;
        vector<double> walker;
        double stima;
        double errore;
        double blk_norm;

        Options function;

};


#endif 