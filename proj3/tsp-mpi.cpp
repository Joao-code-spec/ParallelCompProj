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

bestTaC tspbb(std::vector<std::vector<double>> distances, int nCities, double bestTourCost, int rank, int num_procs){

    std::vector<double> min1 (nCities,INFINITY);
    std::vector<double> min2 (nCities,INFINITY);
    /*used so that we dont have to reshape returnable vector*/
    std::vector<int> rTourFiller (nCities+1,0);
    PriorityQueue<qElement,cmp_op>  queue;
    double d, neiborRet,bestTCResiver;
    bool contains[nCities];
    bool sentFirst;
    int help, rankNext, rankPrev;
    double lowerBound;
    double newBound;
    int token,myColour;
    MPI_Request request, reqForQL[4];
    MPI_Status statsForQL[4];
    MPI_Request reqForReduce;
    MPI_Request reqForBroadc/*[num_procs-1]*/;
    qElement poppedE;
    qElement e;
    vector<int> tour /*= {0}*/;

    char retBuff[270];
    char myBalBuff[280];
    char lnBalBuff[280];
    bool allWhite=false;
    int jkjk=0;
    int step=1;
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
    

    rankNext=(rank + 1)%num_procs;
    rankPrev= rank==0 ? num_procs-1 : rank - 1;
    lowerBound = lb(distances, nCities);
    /*creates the initial state for each queue*/
    help=rank+1;
    while(help<nCities){
        if(distances[0][help]!=INFINITY){
            tour={0,help};
            d = distances[0][help];
            newBound=updateBound(lowerBound, 0, help, min1, min2, distances[0][help]);
            e={tour,d,newBound,2,help};
            queue.push(e);
        }
        help+=num_procs;
    }

    sentFirst=false;
    //token white=0, black=1;
    token=0;
    myColour=0;
    /*qElement e={tour,0,lowerBound,1,0};
    queue.push(e);*/
    bestTaC returnable= {{0},bestTourCost};
    returnable.bt=rTourFiller;
    while(allWhite==false){
        if(queue.empty() != true){
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
                queue.clear();
                //return returnable;
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
                int i=0;
                for(double v : distances[poppedE.currentCity]){
                    if( v != INFINITY && !contains[i]){
                        lowerBound=updateBound(poppedE.bound, poppedE.currentCity, i, min1, min2, distances[poppedE.currentCity][i]);
                        if(lowerBound>bestTourCost){
                            i++;
                            continue;
                        }
                        int newLenght =poppedE.lenght + 1;
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
        /*check neibors every 50 step*/
        if(step%50==0){
            int y=queue.size();
            int x, z;
            qElement eFromPrev={{1}, 2.3, 2.3, 2, 2};
            MPI_Irecv(&x,1,MPI_INT,rankNext,5,MPI_COMM_WORLD,&reqForQL[2]);
            MPI_Irecv(&z,1,MPI_INT,rankPrev,4,MPI_COMM_WORLD,&reqForQL[3]);

            MPI_Isend(&y,1,MPI_INT,rankNext,4,MPI_COMM_WORLD,&reqForQL[0]);
            MPI_Isend(&y,1,MPI_INT,rankPrev,5,MPI_COMM_WORLD,&reqForQL[1]);
            //broadCast version
            MPI_Iallreduce(&bestTourCost,&bestTCResiver,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,&reqForReduce);
            

            
            //Balance Sends one to next if nexts queue is shorter by 20
            /*MPI_Waitall(4,reqForQL,statsForQL);
            if(x<y+20){
                poppedE=queue.pop();
                memcpy(myBalBuff, &poppedE.currentCity, sizeof(int));
                memcpy(myBalBuff+sizeof(int), &poppedE.lenght, sizeof(int));
                memcpy(myBalBuff+2*sizeof(int), &poppedE.bound, sizeof(double));
                memcpy(myBalBuff+2*sizeof(int)+sizeof(double), &poppedE.cost, sizeof(double));
                memcpy(myBalBuff+2*sizeof(int)+2*sizeof(double), poppedE.tour.data(), poppedE.lenght*sizeof(int));
                MPI_Send(myBalBuff,280,MPI_BYTE,rankNext,5,MPI_COMM_WORLD);
                //TODO necesary?
                //turn black if sending to earlier in the ring, only appens at the end of the ring
                if(rankNext==0){
                    myColour=1;
                }
            }

            if(y<=z+20){
                MPI_Recv(lnBalBuff,280,MPI_BYTE,rankPrev,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                memcpy(&eFromPrev.currentCity, lnBalBuff, sizeof(int));
                memcpy(&eFromPrev.lenght, lnBalBuff+sizeof(int), sizeof(int));
                memcpy(&eFromPrev.bound, lnBalBuff+2*sizeof(int), sizeof(double));
                memcpy(&eFromPrev.cost, lnBalBuff+2*sizeof(int)+sizeof(double), sizeof(double));
                eFromPrev.tour.resize(eFromPrev.lenght);
                memcpy(eFromPrev.tour.data(),lnBalBuff+2*sizeof(int)+2*sizeof(double),eFromPrev.lenght*sizeof(int));
                queue.push(eFromPrev);
            }*/

            /*waits for reduce to finish and equalizes all bestTourCosts to the smallest*/
            MPI_Wait(&reqForReduce,MPI_STATUS_IGNORE);
            bestTourCost=bestTCResiver;

            //version send to next
            /*Send first*/
            /*if(rank%2==0){
                MPI_Send((void *)&bestTourCost, 1, MPI_DOUBLE, rankNext, 1, MPI_COMM_WORLD);
                MPI_Recv((void *)&neiborRet, 1, MPI_DOUBLE, rankPrev, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                if(neiborRet<bestTourCost){
                    bestTourCost=neiborRet;
                }
            }*/
            /*Recive first*/
            /*else{
                MPI_Recv((void *)&neiborRet, 1, MPI_DOUBLE, rankPrev, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Send((void *)&bestTourCost, 1, MPI_DOUBLE, rankNext, 1, MPI_COMM_WORLD);
                if(neiborRet<bestTourCost){
                    bestTourCost=neiborRet;
                }
            }*/
            //Termination
	        if (queue.empty()==true) {
                //MPI_Recv(&token, 1, MPI_INT, rankPrev, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int flag=0;
                if (rank == 0) {
                    if(!sentFirst){
                        /*send first token*/
                        MPI_Isend(&token, 1, MPI_INT, rankNext, 2, MPI_COMM_WORLD,&request);
                        sentFirst=true;
                        MPI_Request_free(&request);
                    }
                    else{
                        MPI_Iprobe(rankPrev,2,MPI_COMM_WORLD,&flag,MPI_STATUS_IGNORE);
                        if(flag){
                            MPI_Recv(&token,1,MPI_INT,rankPrev,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                            if(token==0){
                                /*all white terminate process*/
                                allWhite=true;
                            }
                            else{
                                /*sets token to 0 white and restarst cicle*/
                                MPI_Isend(&token, 1, MPI_INT, rankNext, 2, MPI_COMM_WORLD,&request);
                                token=0;
                                MPI_Request_free(&request);
                            }
                        }

                    }
                    
                }
                else {
                    // increment the token and pass it to the next process
                    MPI_Iprobe(rankPrev,2,MPI_COMM_WORLD,&flag,MPI_STATUS_IGNORE);
                    if(flag){
                        MPI_Recv(&token,1,MPI_INT,rankPrev,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        /*black when send to prev node and white when send token to next*/
                        if(myColour==1){
                            token=1;
                        }
                        MPI_Isend(&token, 1, MPI_INT, rankNext, 2, MPI_COMM_WORLD,&request);
                        myColour=0;
                        MPI_Request_free(&request);
                    }
                    //token++;
                    //MPI_Send(&token, 1, MPI_INT, rankNext, 0, MPI_COMM_WORLD);
                }
            }
            MPI_Bcast(&allWhite,1,MPI_CXX_BOOL,0,MPI_COMM_WORLD);
        }
        step++;
    }
    /*makes shore returnable of 0 is the smalest*/
    if(rank==0){
        for(int P=1;P<num_procs;P++){
            MPI_Recv(retBuff,sizeof(double)+((nCities+1)*sizeof(int)),MPI_BYTE,P,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            memcpy(&neiborRet,retBuff,sizeof(double));
            if(neiborRet<returnable.btCost){
                memcpy(returnable.bt.data(),retBuff+sizeof(double),(nCities+1)*sizeof(int));
                returnable.btCost=neiborRet;
            }
        }
    }
    //children send to root
    else{
        memcpy(retBuff,&returnable.btCost,sizeof(double));
        memcpy(retBuff+sizeof(double),returnable.bt.data(),(nCities+1)*sizeof(int));
        MPI_Send(retBuff,sizeof(double)+((nCities+1)*sizeof(int)),MPI_BYTE,0,3,MPI_COMM_WORLD);
    }
    return returnable;
}
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
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

    t=tspbb(roadMatrix,totalCitys,maxVal,rank,num_procs);

    exec_time += omp_get_wtime();
    if(rank==0){
        fprintf(stderr, "%.1fs\n", exec_time);

        if(t.btCost>=maxVal){
            printf("NO SOLUTION\n");
            //return 0;
        }
        else{
            printf("%.1f\n",t.btCost);
            for(int iiii : t.bt){
                printf("%d ",iiii);
            }
            printf("\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
