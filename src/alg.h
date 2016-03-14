#include <string>
#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib> 
#include <fstream>
using namespace std;

// 
//  Change any of these parameters to match your needs 
//
# define POPSIZE 1000//8000
# define MAXGENS 500//1000
# define NVARS 3
# define PXOVER 0.8
# define PMUTATION 0.15

//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
struct genotype
{
	double gene[NVARS];
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	double rfitness;
	double cfitness;
};

struct genotype population[POPSIZE + 1];
struct genotype newpopulation[POPSIZE + 1];

//CPU functions
int main(int argc, char **argv);
void crossover(int &seed);
void elitist();
void evaluate();
int i4_uniform_ab(int a, int b, int &seed);
void initialize(int &seed);
void keep_the_best();
void mutate(int &seed);
double r8_uniform_ab(double a, double b, int &seed);
void report(int generation);
void selector(int &seed);
void timestamp();
void Xover(int one, int two, int &seed);

