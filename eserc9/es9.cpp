#include "eserc9classes.h"
#include "eserc9functions.h"

using namespace std;

int main(int argc, char* argv[]){

    Random rnd;
    initRand(rnd);

    Map m(rnd, 0, 34, 1);

    ofstream DrawCities("Square_map.dat");
    for(int i = 0; i < m.GetNcities(); i++){
        DrawCities << setw(18) << m.GetMap()[i][0] << setw(18) << m.GetMap()[i][1] << endl;
    }
    DrawCities.close();

    int nindividuals = 500;
    int ngenerations = 1500;

    Genetic gen(rnd, m, nindividuals);

    ofstream AverageLengths("AverageLengths.dat");

    for(int i = 0; i < ngenerations; i++){

        gen.Evolve();
        gen.Mutate();

        AverageLengths << setw(10) << i << setw(10) << gen.GetMeanLength() << endl;
    }


    AverageLengths.close();

    ofstream DrawCitiesBest("Cities_square_L1.dat");
    DrawCitiesBest << gen.GetPop()[0].GetLength() << endl;

    for(int i = 0; i < m.GetNcities(); i++){
        DrawCitiesBest << setw(18) << gen.GetPop()[0].GetPath()[i] << endl;
    }
    DrawCitiesBest.close();


//###############################Debug distanze
    // ofstream Distances("circle_distances.dat");
    // for(int i=0; i < m.GetNcities(); i++){

    //     for(int j=0; j < m.GetNcities(); j++){

    //         Distances << setw(15) << m.GetDistances()[i][j];
    //     }
    //     Distances << endl;
    // }
    // Distances.close();
//#####################################

//###############################Debug operatori mutazione
    // Path path(rnd, m);

    // for(int i=0; i < path.GetNcities(); i++){
    //     cout << setw(4) << path.GetPath()[i];
    // }
    // cout << endl << "Lunghezza iniziale: " << path.GetLength() << endl;

    // path.Invert(rnd);
    // path.CalculateLength(m);
    // for(int i=0; i < path.GetNcities(); i++){
    //     cout << setw(4) << path.GetPath()[i];
    // }
    // cout << endl << "Lunghezza finale: " << path.GetLength() << endl;
//###################################


    return 0;

}
