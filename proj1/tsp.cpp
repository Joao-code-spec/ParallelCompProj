#include <stdio.h>
#include "nqueue/queue.h"
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
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
    lb = std::ceil(lb / 2);
    return lb;
}

double updateBound(std::list<int> tour, double cost, int currentCity, int remainingCities, std::vector<std::vector<double>> distances) {
    double lb = cost;
    // Calculate the sum of the two smallest distances for each unvisited city
    for(int i = 0; i < distances.size(); i++) {
        if(std::find(tour.begin(), tour.end(), i) == tour.end()) {
            double min1 = std::numeric_limits<double>::max();
            double min2 = std::numeric_limits<double>::max();
            for(int j = 0; j < distances.size(); j++) {
                if(i != j && std::find(tour.begin(), tour.end(), j) == tour.end()) {
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
    }

    return lb / 2;
}

bestTaC tspbb(std::vector<std::vector<double>> distances, int nCities, double bestTourCost){
    list<int> tour = {0};
    double lowerBound = lb();
    double d;
    qElement e={tour,0,lowerBound,1,0};
    PriorityQueue<qElement>  queue;
    qElement poppedE;
    bestTaC returnable= {{0},9999999999};
    while(queue.empty() != true){
        poppedE=queue.pop();
        if(poppedE.bound>=bestTourCost){
            //poppedE.tour.push_front(poppedE.currentCity);
            returnable={poppedE.tour, bestTourCost};
            return returnable;
        }
        if(poppedE.lenght==nCities){
            //re-used lowerBound because it is a double this has nothing to do with lowerbound
            lowerBound = poppedE.cost + distances[poppedE.currentCity][0];
            if(lowerBound<bestTourCost){
                tour=poppedE.tour;
                tour.push_front(0);
                returnable.bt=tour;
                returnable.btCost=lowerBound;
            }
        }
        else{
            //TODO find better way of finding if contains, look at sets
            int i=0;
            for(double v : distances[poppedE.currentCity]){
                bool contains=false;
                for(int city : poppedE.tour){
                    if(city==i){
                        contains==true;
                        break;
                    }
                }
                if(0 < v && !contains){
                    lowerBound=updateBound(e.tour, e.cost, e.currentCity, nCities - e.lenght, distances);
                    if(lowerBound>bestTourCost){
                        i++;
                        continue;
                    }
                    int i1 =poppedE.lenght + 1;
                    tour=poppedE.tour;
                    tour.push_front(i);
                    d = poppedE.cost + distances[poppedE.currentCity][i];

                    qElement next = {tour,d,lowerBound,i1,v};
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
    int i1, i2, distance;
    if ((file = fopen(argv[1],"r")) == NULL){
       printf("Error! file doesnt exist \n");
       return 1;
    }
    fscanf(file,"%d %d", &totalCitys, &totalRoads);
    //TODO check if broken
    std::vector<std::vector<double>> roadMatrix(totalCitys, std::vector<double>(totalCitys));
    //double roadMatrix [totalCitys][totalCitys];
    //printf("%d %d\n",totalCitys,totalRoads);
    while(fscanf(file,"%d %d %d", &i1, &i2, &distance) != EOF){
        roadMatrix[i1][i2]=distance;
        roadMatrix[i2][i1]=distance;
        //printf("%d %d %d\n", i1, i2, distance);
    }
    fclose(file);
    maxVal = strtol(argv[2], NULL, 10);;
    //printf("%d\n",maxVal);
    return 0;
}
