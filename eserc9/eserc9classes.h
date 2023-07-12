#ifndef __eserc9classes__
#define __eserc9classes__

#include "random.h"
#include "eserc9functions.h"

using namespace std;

class Map {

    public: 

        Map(){}

        Map(Random& rnd, int config, int ncities, int norm) : 
            m_map(ncities, vector<double> (2, 0)), m_distances(ncities, vector<double> (ncities, 0)){

            m_config = config;                 //config=0 is square, else is circle
            m_ncities = ncities;               //number of cities            
            m_norm = norm;                     //norm = 1 or 2 

            if(m_config==0){     //random cities in a square

                for(int i=0; i < m_ncities ; i++){

                    m_map[i][0] = rnd.Rannyu(-1.,1.);
                    m_map[i][1] = rnd.Rannyu(-1.,1.);
                }
            }
            else{    //cities on a circumference  

                for(int i=0; i < m_ncities; i++){

                    double theta = rnd.Rannyu(0., 2.*M_PI);

                    m_map[i][0] = cos(theta);
                    m_map[i][1] = sin(theta);
                    
                }
            }

            //Calculate distances between cities
            for(int i=0; i < m_ncities; i++){       

                for(int j=0; j < m_ncities; j++){

                    if(norm == 1){

                        m_distances[i][j] = pow( pow(m_map[i][0] - m_map[j][0], 2) + pow(m_map[i][1] - m_map[j][1], 2)  , 0.5);

                    }
                    if(norm == 2){

                        m_distances[i][j] = pow(fabs(m_map[i][0] - m_map[j][0]), 2) + pow(fabs(m_map[i][1] - m_map[j][1]), 1);

                    }

                }
            }

        }

        ~Map(){}

        vector<vector<double>> GetMap(){
            return m_map;
        }
        vector<vector<double>> GetDistances(){
            return m_distances;
        }
        int GetConfig(){
            return m_config;
        }
        int GetNorm(){
            return m_norm;
        }
        int GetNcities(){
            return m_ncities;
        }

    private:
        vector<vector<double>> m_map;
        vector<vector<double>> m_distances;
        int m_config, m_norm;
        int m_ncities;

};

class Path {

    public:
        Path(){}

        Path(Random& rnd, Map& map) : m_path(map.GetNcities(), 0) {

            m_ncities = map.GetNcities();

            for(int i = 0; i < m_ncities; i++){       //Riempi con le città in ordine
                m_path[i] = i;
            }
            for(int i=0; i<50; i++){        //Scegli a caso due posizioni da 1 a 33 e scambia i valori (la posizione 0 non può essere cambiata)

                int k = ceil(rnd.Rannyu(0., (double) (m_ncities - 1)));
                int l = ceil(rnd.Rannyu(0., (double) (m_ncities - 1)));
                swap(m_path[k], m_path[l]); 
            }

            double length = 0.;
            for(int i=0; i<m_ncities-1; i++){   //Calcola la lunghezza del cammino

                length += map.GetDistances()[m_path[i]][m_path[i+1]];
            }
            
            length += map.GetDistances()[0][m_path[m_ncities - 1]];       //Aggiunge il tratto tra l'ultima città e la prima.
            m_length = length;
        }

        ~Path(){}

        bool CheckPath(){

            if(m_path[0] != 0){return false;}

            for(int i = 0; i < m_ncities; i++){

                for(int j = i + 1; j < m_ncities; j++){

                    if(m_path[i] == m_path[j]){return false;}
                }   
            }
            return true;
        }

        void CalculateLength(Map& map){

            double length = 0.;
            for(int i=0; i<m_ncities-1; i++){   //Calcola la lunghezza del cammino

                length += map.GetDistances()[m_path[i]][m_path[i+1]];
            }
            length += map.GetDistances()[0][m_path[m_ncities - 1]];       //Aggiunge il tratto tra l'ultima città e la prima.
            m_length = length;
        }

        bool operator<(const Path& p) const{
            return GetLength() < p.GetLength();
        } 

        void SetPath(vector<int> p){
            for(int i=0; i < m_ncities; i++){
                m_path[i] = p[i];
            }
        }

        void Swap(Random& rnd){  //swap the elements with index k and l.

            int k = 0;
            int l = 0;

            while( k==l ){
                k = floor(rnd.Rannyu(1., m_ncities));
                l = floor(rnd.Rannyu(1., m_ncities));
            }
            swap(m_path[k], m_path[l]);
        }
        void Shift(Random& rnd){   //shift by shift positions, num cities, from beg index included onward.

            int beg = ceil(rnd.Rannyu(0., m_ncities-1));      //shift va scelto in modo che gli elementi traslati non si sovrappongano al gruppo iniziale.
            int num = floor(rnd.Rannyu(1., m_ncities-1));     //Quindi deve essere beg + num + shift = beg + map.GetN() ossia num + shift < map.GetN().
            int shift = floor(rnd.Rannyu(1., m_ncities-num)); 

            vector<int> var;      

            if(shift < num){ //Due casi: se lo shift è minore del numero di città, e allora si sovrappongono, oppure no.

                for(int i=0; i < num; i++){     //get first num elements, get them in var
                    var.push_back(m_path[Pbc(m_path, i + beg)]);}
                
                for(int i=0; i < shift; i++){    //get the rest of the elements
                    var.push_back(m_path[Pbc(m_path, i + beg + num)]);}

                rotate(var.begin(), var.begin() + num, var.end());

                for(int i=0; i < shift; i++){
                    m_path[Pbc(m_path, i + beg)] = var[i];}

                for(int i=0; i < num; i++){
                    m_path[Pbc(m_path, i +  beg + shift)] = var[ i + shift ];}
            }
            else{

                for(int i=0; i < num; i++){
                    var.push_back(m_path[Pbc(m_path, i + beg)]);}
                
                for(int i=0; i < num; i++){
                    var.push_back(m_path[Pbc(m_path, i + beg + shift)]);}

                rotate(var.begin(), var.begin() + num, var.end());

                for(int i=0; i < num; i++){
                    m_path[Pbc(m_path, i + beg)] = var[i];}

                for(int i=0; i < num; i++){
                    m_path[Pbc(m_path, i +  beg + shift)] = var[ i + num ];}
            }
        }
        void Rotate(Random& rnd){   //beg and middle are indeces. cities do not overlap here. 

            int beg = ceil(rnd.Rannyu(0., m_ncities - 1));
            int num = floor(rnd.Rannyu(1., (double) m_ncities/2.));    //num is between 1 and N/2 excluded.
            int middle = floor(rnd.Rannyu(beg + num, beg + m_ncities - num)); //numero massimo è tale che middle + num = beg + n. numero minimo è beg + num, senza pbc.

            vector<int> var;

            for(int i=0; i < num; i++){     //get first num elements, get them in var
                var.push_back(m_path[Pbc(m_path, i + beg)]);}
                
            for(int i=0; i < num; i++){    //get the rest of the elements
                var.push_back(m_path[Pbc(m_path, i + middle)]);}

            rotate(var.begin(), var.begin() + num, var.end());

            for(int i=0; i < num; i++){
                m_path[Pbc(m_path, i + beg)] = var[i];}

            for(int i=0; i < num; i++){
                m_path[Pbc(m_path, i + middle)] = var[ i + num ];}  
        }
        void Invert(Random& rnd){  //inverse the order of m contiguous cities

            int beg = ceil(rnd.Rannyu(0., m_ncities - 1));
            int group = ceil(rnd.Rannyu(1., m_ncities - 1));

            vector<int> var;

            for(int i=0; i < group; i++){     //get first num elements, get them in var
                var.push_back(m_path[Pbc(m_path, i + beg)]);}

            reverse(var.begin(), var.end());

            for(int i=0; i < group; i++){
                m_path[Pbc(m_path, i + beg)] = var[i];}
        }

        
        vector<int> GetPath(){
            return m_path;
        }
        int GetNcities(){
            return m_ncities;
        }
        double GetLength() const{
            return m_length;
        }

    private:
        vector<int> m_path;
        int m_ncities;
        double m_length;

};

class Genetic {

    public: 

        Genetic(){}

        Genetic(Random& r, Map& m, int nindividuals) : rnd{r}, map{m}, population(nindividuals) {

            m_nindividuals = nindividuals;

            for(int i = 0; i < m_nindividuals; i++){
                Path p(rnd, map);
                population[i] = p;
            }
            Sort();
        }

        ~Genetic(){}

        void Sort(){
            for(int i = 0; i < m_nindividuals; i++){
                population[i].CalculateLength(map);
            }
            sort(population.begin(),population.end());
        }

        int Select(){

            double p = 4.0;
            return (int) (m_nindividuals * pow(rnd.Rannyu(), p));
        }

        void Mutate(){

            double r1, r2, r3, r4;

            for( int i = 0; i < m_nindividuals; i++){

                r1 = rnd.Rannyu();
                r2 = rnd.Rannyu();
                r3 = rnd.Rannyu();
                r4 = rnd.Rannyu();

                if( r1 < 0.05 ){
                    population[i].Swap(rnd);
                }
                if( r2 < 0.05 ){
                    population[i].Shift(rnd);
                }
                if( r3 < 0.05 ){
                    population[i].Rotate(rnd);
                }
                if( r4 < 0.05 ){
                    population[i].Invert(rnd);
                }
            }
            Sort();
        }

        void Crossover(int selected_1, int selected_2){

            Path child_1 = population[selected_1];
            Path child_2 = population[selected_2];

            double r = rnd.Rannyu();

            if( r < 0.9 ){         //Crossover probability

                int count1 = 0;
                int count2 = 0;

                int cut = floor(rnd.Rannyu(1., map.GetNcities() - 1));      //random cut index

                vector<int> path_1 = child_1.GetPath();
                vector<int> path_2 = child_2.GetPath();

                path_1.resize(cut);             //Paths stay the same until cut point
                path_2.resize(cut);

                // Crossover

                for(int i = 0; i < map.GetNcities(); i ++){    // Iterate over cities in the parent routes

                    count1 = 0;
                    count2 = 0;
                    // Check if the city from selected2 already exists in child1's route
                    for(int j = 0; j < path_1.size() ; j++){
                        if( population[selected_2].GetPath()[i] == path_1[j] ){
                            count1 ++;
                            break;
                        }
                    }
                    // If city doesn't exist, add it to child1's route
                    if (count1 == 0 ) path_1.push_back(population[selected_2].GetPath()[i]);

                    // Check if the city from selected1 already exists in child2's route
                    for(int j = 0; j < path_2.size(); j++){
                        if( population[selected_1].GetPath()[i] == path_2[j] ){
                            count2 ++;
                            break;
                        }
                    }
                    // If city doesn't exist, add it to child2's route
                    if (count2 == 0 ) path_2.push_back(population[selected_1].GetPath()[i]);
                }

                child_1.SetPath(path_1);
                child_2.SetPath(path_2);

            }
            // Add the modified (or not) individuals to the population
            population.push_back(child_1);
            population.push_back(child_2);
        }

        void Evolve(){

            int selected_1 = 0;
            int selected_2 = 0;

            for(int i = 0; i < m_nindividuals / 2; i++){

                selected_1 = Select();
                selected_2 = Select();

                Crossover(selected_1, selected_2);
            }
            Sort();
            population.resize(m_nindividuals);
        }

        
        vector<Path> GetPop(){
            return population;
        }
        int GetNIndividuals(){
            return m_nindividuals;
        }
        double GetMeanLength(){

            double sum = 0.;
            for(int i = 0; i < m_nindividuals/2; i++){
                sum += population[i].GetLength();
            }
            return sum / (double)(m_nindividuals/2.);
        }
    
    private:
        Random rnd;
        Map map;
        vector<Path> population;
        int m_nindividuals;





};





#endif