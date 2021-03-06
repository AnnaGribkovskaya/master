#include <iostream>
#include "quantumdot.h"

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int NumberOfShells = 4;
    int NumberOfElectrons = 2;
    double HOStrenth = 1;


    // If a first argument is provided, it is the number of shells
    if(numberOfArguments > 1) NumberOfShells = atoi(argumentList[1]);
    // If a second argument is provided, it is the number of electrons
    if(numberOfArguments > 2) NumberOfElectrons = atoi(argumentList[2]);
    // If a third argument is provided, it is the HO strenth /omega
    if(numberOfArguments > 3) HOStrenth = atof(argumentList[3]);

    QuantumDot qdot(NumberOfShells, HOStrenth, NumberOfElectrons);
    //qdot.getQuantumDotStates();
    //qdot.applyHartreeFockMethod();
    qdot.applyCoupledClusterDoubles();
    //qdot.getQuantumDotStatesNumber();
    //qdot.applyCCD_Not_HF();


}
