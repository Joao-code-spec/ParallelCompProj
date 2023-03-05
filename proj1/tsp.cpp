#include <stdio.h>
#include "nqueue/queue.h"
#include "nqueue/queue.hpp"
#include <stdlib.h>
#include <iostream> 
#include <list>
using namespace std;
int maxVal;
typedef struct
{
	list<int> tour;
    double cost;
	double lb;
    int lenght;
    int currentCity;
} qElement;

double lb(){
    //TODO
    return 0;
}
int tspbb(double distances, int nCities, double bestTourCost){
    list<int> tour = {0};
    double lowerBound = lb();
    qElement e={tour,0,lowerBound,1,0};
    PriorityQueue<qElement>  queue;

    qElement poppedQueue;
    while(queue.empty() != true){
        poppedQueue=queue.pop();
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