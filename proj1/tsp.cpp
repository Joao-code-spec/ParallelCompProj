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
	list<int> tour;
    double cost;
	double bound;
    int lenght;
    int currentCity;
} qElement;
struct cmp_op { bool operator()(qElement const&  left,qElement const& right) { return left.bound > right.bound; } };

typedef struct
{
	list<int> bt;
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

double updateBound(double cost, int currentCity, std::vector<double> distanceCurrent,std::vector<double> distanceTo, double dct) {
    double lb = cost;
    double cf,ct,cft;
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

    }
    if(dct>=min2f){
        cf=min2f;
    }
    else{
        cf=min1f;
    }
    if(dct>=min2t){
        ct=min2t;
    }
    else{
        ct=min1t;
    }
    cft=cf +ct;
    cft= cft /2;
    return lb + dct - cft;
}

bestTaC tspbb(std::vector<std::vector<double>> distances, int nCities, double bestTourCost){
    list<int> tour = {0};
    double lowerBound = lb(distances, nCities);
    double d;
    qElement e={tour,0,lowerBound,1,0};
    PriorityQueue<qElement,cmp_op>  queue;
    queue.push(e);
    qElement poppedE;
    bestTaC returnable= {{0},bestTourCost};
    while(queue.empty() != true){
        poppedE=queue.pop();
        //printf("POPPED %d LB %.1f cost %.1f \n",poppedE.currentCity,poppedE.bound,poppedE.cost);
        if(poppedE.bound>=bestTourCost){
            
            //poppedE.tour.push_front(poppedE.currentCity);
            returnable={returnable.bt, bestTourCost};
            return returnable;
        }
        if(poppedE.lenght==nCities){
            //re-used lowerBound because it is a double this has nothing to do with lowerbound
            // if Cost + Distances(Node, 0) < BestT ourCost then
            lowerBound = poppedE.cost + distances[poppedE.currentCity][0];
            if(lowerBound<bestTourCost){
                tour=poppedE.tour;
                tour.push_back(0);
              	returnable.bt=tour;
                returnable.btCost=lowerBound;
                bestTourCost=lowerBound;
            }
        }
        else{
            //TODO find better way of finding if contains, look at sets
            int i=0;
            for(double v : distances[poppedE.currentCity]){
                bool contains=false;
                for(int city : poppedE.tour){
                    if(city==i){
                        contains=true;
                        break;
                    }
                }
                if( v != INFINITY && !contains){
                    lowerBound=updateBound(poppedE.bound, poppedE.currentCity, distances[poppedE.currentCity],distances[i],distances[poppedE.currentCity][i]);
                    if(lowerBound>bestTourCost){
                        i++;
                        continue;
                    }
                    int newLenght =poppedE.lenght + 1;
                    tour=poppedE.tour;
                    tour.push_back(i);
                    d = poppedE.cost + distances[poppedE.currentCity][i];

                    qElement next = {tour,d,lowerBound,newLenght,i};
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
    //exec_time = -omp_get_wtime();

    t=tspbb(roadMatrix,totalCitys,maxVal);

    //exec_time += omp_get_wtime();
    //fprintf(stderr, "%.1fs\n", exec_time);

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
