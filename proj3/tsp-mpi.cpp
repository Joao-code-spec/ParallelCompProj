#include <stdio.h>
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <mpi.h>

using namespace std;
int maxVal;
typedef struct
{
	vector<int> tour;
    double cost;
	double bound;
    int lenght;
    int currentCity;
} qElement;

struct cmp_op { bool operator()(qElement const&  left,qElement const& right) { 
    if(left.bound == right.bound){
        return left.currentCity>right.currentCity;
    }
    return left.bound > right.bound; } };

typedef struct
{
	vector<int> bt;
	double btCost;
} bestTaC;

    
double lb(std::vector<std::vector<double>> distances, int nCities) {
    double lb = 0;
    // calculate the sum of the two smallest distances for each city
    for(int i = 0; i < nCities; i++) {
        double min1 = std::numeric_limits<double>::max();
        double min2 = std::numeric_limits<double>::max();
        for(int j = 0; j < nCities; j++) {
            if(i != j) {
                double dist = distances[i][j];
                if(dist < min1) {
                    min2 = min1;
                    min1 = dist;
                } else if(dist < min2) {
                    min2 = dist;
                }
            }
        }
        lb += min1 + min2;
    }
    // divide the result by 2 and round up to the nearest integer
    lb = lb / 2;
    return lb;
}


double updateBound(double cost, int currentCity,int toCity,vector<double> &min1,vector<double> &min2, double dct) {
    double lb = cost;
    double cf,ct,cft;
    /*
    double min1f = INFINITY;
    double min2f = INFINITY;
    double min1t = INFINITY;
    double min2t = INFINITY;
    // Calculates mins for current city
    for(int akf = 0; akf < (int) distanceCurrent.size();akf++){
        double dist = distanceCurrent[akf];
        if(dist < min1f) {
            min2f = min1f;
            min1f = dist;
        } else if(dist < min2f) {    
            min2f = dist;
        }

    }
    //calculates mins for toCity
    for(int akf = 0; akf < (int) distanceTo.size();akf++){
        double dist = distanceTo[akf];
        if(dist < min1t) {
            min2t = min1t;
            min1t = dist;
        } else if(dist < min2t) {    
            min2t = dist;
        }

    }*/
    if(dct>=min2[currentCity]){
        cf=min2[currentCity];
    }
    else{
        cf=min1[currentCity];
    }
    if(dct>=min2[toCity]){
        ct=min2[toCity];
    }
    else{
        ct=min1[toCity];
    }
    cft=cf +ct;
    cft= cft /2;
    return lb + dct - cft;
}

bestTaC tspbb(std::vector<std::vector<double>> distances, int nCities, double bestTourCost){

    std::vector<double> min1 (nCities,INFINITY);
    std::vector<double> min2 (nCities,INFINITY);
    double d;
    bool contains[nCities];
    int help;
    vector<int> tour = {0};
    int jkjk=0;
    for(std::vector<double> cdd : distances){
        for(int akf = 0; akf < (int) cdd.size();akf++){
            //double dist = cdd[akf];
            if(cdd[akf] < min1[jkjk]) {
                min2[jkjk] = min1[jkjk];
                min1[jkjk] = cdd[akf];
            } else if(cdd[akf] < min2[jkjk]) {    
                min2[jkjk] = cdd[akf];
            }

        }
        jkjk++;
    }
    
    double lowerBound = lb(distances, nCities);
    qElement e={tour,0,lowerBound,1,0};
    PriorityQueue<qElement,cmp_op>  queue;
    queue.push(e);
    qElement poppedE;
    bestTaC returnable= {{0},bestTourCost};

    while(queue.empty() != true){
        poppedE=queue.pop();
        help=0;
        while(help<nCities){
            contains[help]=false;
            help++;
        }
        for(int c : poppedE.tour){
            contains[c]=true;
        }
        //printf("POPPED %d LB %.1f cost %.1f \n",poppedE.currentCity,poppedE.bound,poppedE.cost);
        if(poppedE.bound>=bestTourCost){
            
            //poppedE.tour.push_front(poppedE.currentCity);
            //returnable={returnable.bt, bestTourCost};
            return returnable;
        }
        if(poppedE.lenght==nCities){
            //re-used lowerBound because it is a double this has nothing to do with lowerbound
            // if Cost + Distances(Node, 0) < BestT ourCost then
            lowerBound = poppedE.cost + distances[poppedE.currentCity][0];
            if(lowerBound<bestTourCost){
                //tour=poppedE.tour;
                //tour.push_back(0);
              	returnable.bt=poppedE.tour;
                returnable.bt.push_back(0);
                returnable.btCost=lowerBound;
                bestTourCost=lowerBound;
            }
        }
        else{
            //TODO find better way of finding if contains, look at sets
            int i=0;
            for(double v : distances[poppedE.currentCity]){
                if( v != INFINITY && !contains[i]){
                    lowerBound=updateBound(poppedE.bound, poppedE.currentCity, i, min1, min2, distances[poppedE.currentCity][i]);
                    if(lowerBound>bestTourCost){
                        i++;
                        continue;
                    }
                    int newLenght =poppedE.lenght + 1;
                    //tour=poppedE.tour;
                    //tour.push_back(i);
                    d = poppedE.cost + distances[poppedE.currentCity][i];

                    qElement next = {poppedE.tour,d,lowerBound,newLenght,i};
                    next.tour.push_back(i);
                    //printf("pushed %d LB %.2f cost %.2f \n",next.currentCity,next.bound,next.cost);
                    queue.push(next);
                }
                i++;
            }
        }
    }
    return returnable;
}
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    FILE * file;
    int totalCitys;
    int totalRoads;
    int i1, i2;
    double exec_time, distance;
    bestTaC t;
    if ((file = fopen(argv[1],"r")) == NULL){
       printf("Error! file doesnt exist \n");
       return 0;
    }
    fscanf(file,"%d %d", &totalCitys, &totalRoads);
    //TODO check if broken
    std::vector<std::vector<double>> roadMatrix(totalCitys, std::vector<double>(totalCitys,INFINITY));
    //double roadMatrix [totalCitys][totalCitys];
    //printf("%d %d\n",totalCitys,totalRoads);
    while(fscanf(file,"%d %d %lf", &i1, &i2, &distance) != EOF){
        roadMatrix[i1][i2]=distance;
        roadMatrix[i2][i1]=distance;
        //printf("%d %d %lf\n", i1, i2, distance);
    }
    fclose(file);
    maxVal = strtol(argv[2], NULL, 10);;
    //printf("%d\n",maxVal);
    exec_time = -omp_get_wtime();

    t=tspbb(roadMatrix,totalCitys,maxVal);

    exec_time += omp_get_wtime();
    if(rank==0)
        fprintf(stderr, "%.1fs\n", exec_time);

    if(t.btCost>=maxVal){
        std::cout << "NO SOLUTION\n" << std::endl;
        return 0;
    }
    printf("%.1f\n",t.btCost);
    for(int iiii : t.bt){
        printf("%d ",iiii);
    }
    printf("\n");
    
    MPI_Finalize();
    return 0;
}
