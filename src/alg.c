#include "alg.h"

int main(int argc, char **argv) {
	int generation;
	int i;
	int seed;

	//CPU
	clock_t tStart = clock();

	int rank, size;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   printf ("process %d, size %d\n", rank, size);

	if(rank==0) {
	//проверка на количество неизвестных
	if (NVARS < 2)
	{
		cout << "\n";
		cout << "  —крещивание невозможно\n";
	}

	cout << "POPSIZE=" << POPSIZE << "\n";
	}

	//seed = 123456789;
	seed = abs(12345 * abs(rand()) + rank*abs(rand()));
	//cout << "seed=" << seed << "\n";

	initialize(seed);

	evaluate();

	keep_the_best();

	for (generation = 0; generation < MAXGENS; generation++)
	{
		selector(seed);
		crossover(seed);
		mutate(seed);
		//report(generation);
		evaluate();
		elitist();
	}

	cout << "Rank: "<< rank;
	cout << "  Best member after " << MAXGENS << " generations: ";
	//cout << "\n";

	for (i = 0; i < NVARS; i++)
	{
		cout << " var(" << i << ") = " << population[POPSIZE].gene[i];
	}

	//cout << "\n";
	cout << " Best fitness = " << population[POPSIZE].fitness << " ";
	//
	//  «авершение
	//

	//timestamp();
    cout << "CPU time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << "\n";
//	printf("Rank %d CPU time: %.2fs\n", rank, (double)(clock() - tStart) / CLOCKS_PER_SEC);

	MPI_Finalize();
	return 0;
}

void crossover(int &seed)

//****************************************************************************80
// —крещивание. ¬ыбираютс¤ 2 родител¤ дл¤ одного срещивани¤
//
//   int FIRST счетчик числа выбранных хромосом
//
//    int &SEED стартовое значение генератора случайных чисел
//
{
	const double a = 0.0;
	const double b = 1.0;
	int mem;
	int one;
	int first = 0;
	double x;


	for (mem = 0; mem < POPSIZE; ++mem)
	{
		//debug
		//cout << "mem=" << mem << " POPSIZE=" << POPSIZE << "\n";

		x = r8_uniform_ab(a, b, seed);

		if (x < PXOVER)
		{
			++first;

			if (first % 2 == 0)
			{
				Xover(one, mem, seed);
			}
			else
			{
				one = mem;
			}

		}
	}
	return;
}

void elitist()

//
//    ELITIST сохран¤ет лучшие хромосомы предыдущей попул¤ции
//
//	Ћучшие хромосомы предыдущей генерации сохран¤ютс¤ последними в массиве. 
//	≈сли лучшие хромосомы текущей генерации хуже, чем лучшие в предыдущей, то перва¤ будет заменена худшей в текущей попул¤ции
//
{
	int i;
	double best;
	int best_mem;
	double worst;
	int worst_mem;

	best = population[0].fitness;
	worst = population[0].fitness;


	for (i = 0; i < POPSIZE - 1; ++i)
	{
		if (population[i + 1].fitness < population[i].fitness)
		{

			if (best <= population[i].fitness)
			{
				best = population[i].fitness;
				best_mem = i;
			}

			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}

		}
		else
		{

			if (population[i].fitness <= worst)
			{
				worst = population[i].fitness;
				worst_mem = i;
			}

			if (best <= population[i + 1].fitness)
			{
				best = population[i + 1].fitness;
				best_mem = i + 1;
			}

		}

	}

//	≈сли лучшие индивиды из новой попул¤ции лучше чем лучшие из предыдущей попул¤ции,
//	то копируютс¤ лучшие из новой попул¤ции или замен¤ютс¤ худшие из текущей попул¤ции на лучший из предыдущей попул¤ции

	if (population[POPSIZE].fitness <= best)
	{
		for (i = 0; i < NVARS; i++)
		{
			population[POPSIZE].gene[i] = population[best_mem].gene[i];
		}
		population[POPSIZE].fitness = population[best_mem].fitness;
	}
	else
	{
		for (i = 0; i < NVARS; i++)
		{
			population[worst_mem].gene[i] = population[POPSIZE].gene[i];
		}
		population[worst_mem].fitness = population[POPSIZE].fitness;
	}

	return;
}


void evaluate()

//    EVALUATE функци¤ эволюции
//     x[1]^2-x[1]*x[2]+x[3]
{
	int member;
	int i;
	double x[NVARS + 1];

	for (member = 0; member < POPSIZE; member++)
	{
		for (i = 0; i < NVARS; i++)
		{
			x[i + 1] = population[member].gene[i];
		}
		population[member].fitness = (x[1] * x[1]) - (x[1] * x[2]) + x[3];
	}
	return;
}

int i4_uniform_ab(int a, int b, int &seed)

//    I4_UNIFORM_AB возвращает I4 в диапазоне от A до B
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	r = (float)(seed)* 4.656612875E-10;

	r = (1.0 - r) * ((float)a - 0.5)
		+ r   * ((float)b + 0.5);

	value = round(r);

	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}


void initialize(int &seed)

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE инициализаци¤ хромосом из файла
//		‘ормат файла:
//	нижн¤¤_граница   верхн¤¤_граница
{
	int i;
	int j;

	int ii=0;

	double input[3][2] =  {{0.0, 5.0}, {0.0, 5.0}, {-2.0, 2.0}};

	for (i = 0; i < NVARS; i++)
	{
		for (j = 0; j < POPSIZE; j++)
		{
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].lower[i] = input[ii][0];
			population[j].upper[i] = input[ii][1];
			population[j].gene[i] = r8_uniform_ab(input[ii][0], input[ii][1], seed);
		}
	ii++;
	}

	return;
}

void keep_the_best()

//    KEEP_THE_BEST сохран¤ет текущие хромосомы попул¤ции

{
	int cur_best;
	int mem;
	int i;

	cur_best = 0;

	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[POPSIZE].fitness < population[mem].fitness)
		{
			cur_best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	// 
	//   огда лучшие из попул¤ции найдены - копируем их
	//
	for (i = 0; i < NVARS; i++)
	{
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
	}

	return;
}


void mutate(int &seed)

//    MUTATE случайна¤ мутаци¤

{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	double lbound;
	double ubound;
	double x;

	for (i = 0; i < POPSIZE; i++)
	{
		for (j = 0; j < NVARS; j++)
		{
			x = r8_uniform_ab(a, b, seed);
			if (x < PMUTATION)
			{
				lbound = population[i].lower[j];
				ubound = population[i].upper[j];
				population[i].gene[j] = r8_uniform_ab(lbound, ubound, seed);
			}
		}
	}

	return;
}

double r8_uniform_ab(double a, double b, int &seed)

//    R8_UNIFORM_AB возвращает псевдослучайное R8 в диапазоне от a до b
{
	int i4_huge = 2147483647;
	int k;
	double value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	value = (double)(seed)* 4.656612875E-10;

	value = a + (b - a) * value;

	return value;
}

void report(int generation)

//	выводит значени¤ в процессе работы
{
	double avg;
	double best_val;
	int i;
	double square_sum;
	double stddev;
	double sum;
	double sum_square;

	if (generation == 0)
	{
		cout << "\n";
		cout << "  Generation       Best            Average       Standard \n";
		cout << "  number           value           fitness       deviation \n";
		cout << "\n";
	}

	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < POPSIZE; i++)
	{
		sum = sum + population[i].fitness;
		sum_square = sum_square + population[i].fitness * population[i].fitness;
	}

	avg = sum / (double)POPSIZE;
	square_sum = avg * avg * POPSIZE;
	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
	best_val = population[POPSIZE].fitness;

	cout << "  " << setw(8) << generation
		<< "  " << setw(14) << best_val
		<< "  " << setw(14) << avg
		<< "  " << setw(14) << stddev << "\n";

	return;
}

void selector(int &seed)

//    SELECTOR функци¤ селекции

{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int mem;
	double p;
	double sum;
	//
	//  ¬ычисление общей живучести попул¤ции
	//
	sum = 0.0;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		sum = sum + population[mem].fitness;
	}
	//
	//  –асчет относительной живучести каждого члена попул¤ции
	//
	for (mem = 0; mem < POPSIZE; mem++)
	{
		population[mem].rfitness = population[mem].fitness / sum;
	}
	// 
	//  –асчет совокупной относительной живучести
	//
	population[0].cfitness = population[0].rfitness;
	for (mem = 1; mem < POPSIZE; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness +
			population[mem].rfitness;
	}
	// 
	//  ¬ыбор наиболее приспособленных
	//
	for (i = 0; i < POPSIZE; i++)
	{
		p = r8_uniform_ab(a, b, seed);
		if (p < population[0].cfitness)
		{
			newpopulation[i] = population[0];
		}
		else
		{
			for (j = 0; j < POPSIZE; j++)
			{
				if (population[j].cfitness <= p && p < population[j + 1].cfitness)
				{
					newpopulation[i] = population[j + 1];
				}
			}
		}
	}
	// 
	//  ѕерезаписать старую попул¤цию новой
	//
	for (i = 0; i < POPSIZE; i++)
	{
		population[i] = newpopulation[i];
	}

	return;
}


void Xover(int one, int two, int &seed)

//    XOVER скрещивание 

{
	int i;
	int point;
	double t;
	// 
	//  ¬ыбор точки скрещивани¤ случайным образом
	//
	point = i4_uniform_ab(0, NVARS - 1, seed);
	//cout << "seed=" << seed << "\n";
	//cout << "point=" << point << " one=" << one << " two=" << two << " seed=" << seed << "\n";
	//
	//	ѕомен¤ть местами хромосомы от 0 до POINT-1
	//
	for (i = 0; i < point; i++)
	{
		t = population[one].gene[i];
		//debug
		//cout << "i=" << i << " t=" << t << "\n";
		//cout << "old values:\n" << "population[" << one << "].gene[" << i << "] = " << population[one].gene[i] << "\n";
		//cout << "population[" << two << "].gene[" << i << "] = " << population[two].gene[i] << "\n";

		population[one].gene[i] = population[two].gene[i];
		population[two].gene[i] = t;
		//debug
		//cout << "new values:" << "\n";
		//cout << "population[" << one << "].gene[" << i << "] = " << population[one].gene[i] << "\n";
		//cout << "population[" << two << "].gene[" << i << "] = " << population[two].gene[i] << "\n";

	}

	return;
}
