#include <stdio.h>
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <omp.h>

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
//TODO index of city to break draws as stated in 3.3
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
    bool contains[nCities];
    bool qConfirmedEmpty=false;
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
    //sets one priority queue for each process         numOfthreads
    std::vector<PriorityQueue<qElement,cmp_op>>  queues(1);
    queues[0].push(e);
    PriorityQueue<qElement,cmp_op>  masterQueue;
    qElement poppedE;
    bestTaC returnable= {{0},bestTourCost};


    #pragma omp parallel private(poppedE,lowerBound,contains) num_threads(1)
    {
        /*make step chared by all treads*/
        int step=0;
        int id = omp_get_thread_num();
        int nOfThreads = omp_get_num_threads();
        while(!qConfirmedEmpty){
            if(!queues[id].empty()){
                poppedE=queues[id].pop();
                int help=0;
                while(help<nCities){
                    contains[help]=false;
                    help++;
                }
                for(int c : poppedE.tour){
                    contains[c]=true;
                }
                if(poppedE.bound>=bestTourCost){
                    queues[id].clear();
                }
                if(poppedE.lenght==nCities){
                    //re-used lowerBound because it is a double this has nothing to do with lowerbound
                    // if Cost + Distances(Node, 0) < BestT ourCost then
                    lowerBound = poppedE.cost + distances[poppedE.currentCity][0];
                    #pragma omp critical
                    {
                        if(lowerBound<bestTourCost){
                            returnable.bt=poppedE.tour;
                            returnable.bt.push_back(0);
                            returnable.btCost=lowerBound;
                            bestTourCost=lowerBound;
                        }
                    }
                }
                else{
                    int i=0;
                    for(double v : distances[poppedE.currentCity]){
                        if( v != INFINITY && !contains[i]){
                            lowerBound=updateBound(poppedE.bound, poppedE.currentCity, i, min1, min2, distances[poppedE.currentCity][i]);
                            if(lowerBound>bestTourCost){
                                i++;
                                continue;
                            }
                            int newLenght =poppedE.lenght + 1;
                            double d = poppedE.cost + distances[poppedE.currentCity][i];

                            qElement next = {poppedE.tour,d,lowerBound,newLenght,i};
                            next.tour.push_back(i);
                            queues[id].push(next);
                        }
                        i++;
                    }
                }
            }
            if(step%100 == 0){
                /*waits to merge*/
                #pragma omp barrier
                #pragma omp single
                {
                    /*all empty ?*/
                    bool allEmpty=true;
                    for(PriorityQueue<qElement,cmp_op> q : queues){
                        if(!q.empty()){
                            allEmpty=false;
                        }
                    }
                    if(allEmpty){
                        qConfirmedEmpty=true;
                    }
                    /*merge*/
                    for(PriorityQueue<qElement,cmp_op> q : queues){
                        while(q.empty()!=true){
                            masterQueue.push(q.pop());
                        }
                    }
                    int zx=0;
                    while(!masterQueue.empty()){
                        queues[zx%nOfThreads].push(masterQueue.pop());
                    }
                }
                #pragma omp barrier
            }
        }
        step++;
    }
    return returnable;
}
int main(int argc, char *argv[]){
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
    return 1;
}
