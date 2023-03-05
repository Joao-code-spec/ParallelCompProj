#include <stdio.h>
#include "nqueue/queue.h"
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>
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

double lb(){
    //TODO
    return 0;
}
double updateBound(){
    //TODO
    return 0;
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
                    lowerBound=updateBound();
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