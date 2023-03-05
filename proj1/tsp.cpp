#include <stdio.h>
#include "nqueue/queue.h"
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream> 
#include <list>
using namespace std;
int maxVal;
double lb(){
    //TODO
    return 0;
}
int tspbb(double distances, int nCities, double bestTourCost){
    list<int> tour = {0};
    double lowerBound = lb();
    priority_queue_t * queue = queue_create(tour,0,lowerBound,1,0);
    while(queue != NULL){
        
    }
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
    double roadMatrix [totalCitys][totalCitys];
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