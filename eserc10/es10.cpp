#include "eserc10classes.h"
#include "eserc10functions.h"

using namespace std;

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    //Mappa 
    ifstream read("American_capitals.in");

    vector< vector<double> > citiesIn (50, vector<double> (2,0.));

    for(int i = 0; i < 50; i++){

        read >> citiesIn[i][0] >> citiesIn[i][1];
    }
    read.close();

    Map map(50,1);
    map.SetMap(citiesIn);
    map.CalculateDistances();

    ofstream DrawCities("Map.out");             
    for(int i = 0; i < map.GetNcities(); i++){
        DrawCities << setw(18) << map.GetMap()[i][0] << setw(18) << map.GetMap()[i][1] << endl;
    }
    DrawCities.close();         

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Inizializzazione generatore di numeri casuali. Uno condiviso, gli altri per singolo processo.
    Random shared_rnd;
    initRand(shared_rnd,0);

    Random rnd;
    initRand(rnd,rank);

    //Parametri simulazione
    int nindividuals = 150;
    int ngenerations = 1000;
    int migrate = 1; 
    int N_migr = 10;

    vector<int> sending_order;       //Vettore di indici per lo scambio
    for(int i=0; i < size; i++){
        sending_order.push_back(i);
    }
    for(int i=0; i < 10; i++){
        Shuffle(shared_rnd,sending_order);
    }
 
    int current_rank, send_to, received_from;
    Path best;
    vector<int> best_path;

    Genetic gen(rnd, map, nindividuals);

    
    if(migrate == 0){   //Senza migrazione

        ofstream AverageLengths("AverageLengths_nomig_" + to_string(rank) +  ".dat");

        for(int i = 0; i < ngenerations; i++){

            gen.Evolve();
            gen.Mutate();

            AverageLengths << setw(10) << i << setw(10) << gen.GetMeanLength() << endl;
                
        }
        AverageLengths.close();

    }else{     //With migration

        ofstream AverageLengths("AverageLengths_mig_" + to_string(rank) +  ".dat");

        for(int i = 0; i < ngenerations; i++){

            gen.Evolve();
            gen.Mutate();

            AverageLengths << setw(10) << i << setw(10) << gen.GetMeanLength() << endl;

            if( i!=0 and i%N_migr==0 ){   //exchange individuals every N_migr generations

                for(int j = 0; j < size; j++){        

                    current_rank = sending_order[Pbcnew(sending_order, j)];  
                    send_to = sending_order[Pbcnew(sending_order, j+1)];  
                    received_from = sending_order[Pbcnew(sending_order, j-1)];
                
                    if(rank == current_rank){
                        
                        best = gen.GetPop()[0];
                        best_path = best.GetPath();
                        MPI_Request request1;
                        MPI_Request request2;
                
                        //MPI_Send(best_path.data(), best.GetNcities(), MPI_INT, send_to, 0, MPI_COMM_WORLD);
                        //MPI_Recv(best_path.data(), best.GetNcities(), MPI_INT, received_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        MPI_Isend(best_path.data(), best.GetNcities(), MPI_INT, send_to, 0, MPI_COMM_WORLD, &request1);
                        MPI_Wait(& request1, MPI_STATUS_IGNORE);

                        MPI_Irecv(best_path.data(), best.GetNcities(), MPI_INT, received_from, 0, MPI_COMM_WORLD, & request2);
                        MPI_Wait(& request2, MPI_STATUS_IGNORE);
                  
                        best.SetPath(best_path);
                        gen.SetIndividual(best, gen.GetNIndividuals() - 1);        //sostituisci individuo peggiore con nuovo individuo
        
                    }

                }

                for(int j = 0; j < 10; j++){
                    Shuffle(shared_rnd, sending_order);
                }

                gen.Sort();
            }
            
        }

        AverageLengths.close();
    
    }

    best = gen.GetPop()[0];

    // Declare send and receive buffers for the gathering operation
    vector<int> send = best.GetPath();
    vector<int> receive(size * best.GetNcities());

    // Use MPI_Gather() to collect the information of the best individuals from all nodes
    MPI_Gather(send.data(), best.GetNcities(), MPI_INT, receive.data(), best.GetNcities(), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {    // Root process

        // Population to store the best individuals from each node
        Genetic population_best(shared_rnd, map, size);

        for (int i = 0; i < size; i++) {

            best.SetPath(vector<int>(receive.begin() + (i * best.GetNcities()), receive.begin() + ((i + 1) * best.GetNcities())));
            population_best.SetIndividual(best, i);
        }

        population_best.Sort();
        
        ofstream DrawCitiesBest("Cities.dat");              //Print the best individual
        DrawCitiesBest << gen.GetPop()[0].GetLength() << endl;

        for(int i = 0; i < map.GetNcities(); i++){
            DrawCitiesBest << setw(18) << population_best.GetPop()[0].GetPath()[i] << endl;
        }
        DrawCitiesBest.close();

    }

    MPI_Finalize();


    return 0;

}
