BIOS - Bio-Inspired Optimization System
Version: 4.0.0
-----------------------------------------------------------------------
This program obtains optimum solutions for engineering problems using
a Genetic Algorithm with a hybrid computational parallelization scheme. 
Below is a set of instructions showing how to compile and run the 
program, as well as guidelines on how to build the required data file 
and how to read the obtained results.
-----------------------------------------------------------------------

-----------------------------------------------------------------------
COMPILING
-----------------------------------------------------------------------

To compile the program, one may use the provided makefile (for unix-
based systems only):

Clean old binary and object files:

$ cd obj
$ rm *.o
$ cd ..
$ rm ./bios

Make a clean binary based on the most recent code:

$ make

If gcc is not installed, version 4.7 or later must be first obtained.
In Ubuntu, this can be done by running:

$ sudo apt-get update
$ sudo apt-get install build-essential

Alternatively, one may compile the program in Windows based systems
by including the .cpp and .h files in the folders ctrl, ind, main,
mlib, pen and prob in a project.

The makefile can be modified in order to control whether the two
types of computational parallelization will be used:

To turn on MPI parallelization (turned off by default):
- Remove comment from line 25: CC = mpiCC
- Comment line 26: CC = g++
- Remove comment from line 33: DEFINE = -D_UNIX_ -D_MPI_
- Comment line 34: DEFINE = -D_UNIX_

To turn off OpenMP parallelization (turned on by default):
- Remove the flag '-fopenmp' from line 37.
- Remove the flag '-fopenmp' from line 92.

One can also set the number of cores used in the program execution
by running:

$ export OMP_NUM_THREADS=X, where X is the desired number of cores.


-----------------------------------------------------------------------
RUNNING
-----------------------------------------------------------------------

To run the program, the user may provide the input data file as an argu-
ment:

$ ./bios data_file

or simply run:

$ ./bios

and wait for the program to ask for the file. The program will then be
executed. Is is important to remember that the data file must have the
.opt extension.

To run the program using the Island Model parallelization (based on MPI),
a host file with the name 'mpi_hosts' has to be created. The file used to
run BIOS in the SGI cluster in LMCV (Laboratório de Mecânica Computacional
e Visualização) is provided for reference. Put the hostname of each
computing node in a separate line. It is important to assure ssh access
between all nodes and also between each node and the head node which is
used to run the program. Such direct access may be achieved by generating
ssh authentication keys between all machines.

The file runBIOS.sh must also be edited, specifically the field BIOSFOLDER,
which must contain the path to where the bios binary is located, and the
field FILE, which must contain the path to the data file that will be
executed.

The program can then be executed by running:

$ sh runBIOS.sh

It is important to note that the use of the MPI parallelization scheme
can lead to significant gains in both processing speed and result
reliability. However, the program is also suited to run on common multi-
core and single-core machines.

-----------------------------------------------------------------------
BUILDING THE INPUT FILE
-----------------------------------------------------------------------

The program uses an input file format based on labels, which are marked
by the symbol %. The program then searches every % and identifies the
label. The relevant data must then be located under each label. Each
one of the labels is explained below, where the following conventions 
are adopted:

[x] - x has to be an integer (e.g. 2)
<x> - x has to be real (e.g. 2.25)
'x' - x has to be a char string (e.g. TsaiWu)

%CONSTRAINT.TOLERANCE
<tol>

Sets the desired tolerance when the expression g(x) <= 0 is tested to
determine if a constraint is violated, active or inactive.

%CROSSOVER.RANGE
<MinCross> <MaxCross>

(Values between 0 and 1)
Sets the range for the crossover rate when using the Island Model. Each
island (node) will then have a random rate value between MinCross and
MaxCross.

%CROSSOVER.RATE
<CrossRate>

(Value between 0 and 1)
Sets the crossover rate. The number of individuals selected for crossover
will then be parsize = CrossRate*PopSize, rounded down, and the number of
individuals generated in the crossover process will be sonsize = parsize -
parsize%2.

%INDIVIDUAL.TYPE
'type'

-
IntegerVector
DoubleVector
BinaryVector
IntegerMatrix
-

Chooses the type of individual.

%MAXIMUM.GENERATIONS
[MaxGen]

(Value greater than 0)
Specifies the number of generations for which the genetic algorithm will
run. In general, raising the number of generations will increase the
algorithm reliability.

%MIGRATION.GENERATION.GAP
[MigrGap]

(Value greater than 0)
Sets the generation interval after which the migration process will occur.
When migration occurs, the best individuals of each island are transferred
to all other islands.

%MIGRATION.INDIVIDUAL.NUMBER
[MigrNum]

(Value greater than 0)
Sets the number of individuals of each island which are transferred to
the other islands when migration occurs.

%MUTATION.PROBABILITY
<MutProb>

(Value between 0 and 1)
Sets the desired probability for the Mutation genetic operator. This
operator changes the value of one gene to a random value inside the list
of possible values.

%MUTATION.RANGE
<MinMut> <MaxMut>

(Values between 0 and 1) 
Sets the range for the mutation probability when using the Island Model.
Each island (node) will then have a random probability value between
MinMut and MaxMut.

%OPTIMIZATION.ALGORITHM
'alg'

-
StdGA
-

Sets the desired algorithm which will run the problem.


%OPTIMIZATION.NUMBER
[OptNum]

(Value greater than 0)
Sets the number of times the optimization process will be executed.
In the end of all optimizations, the best obtained individual will be
printed, as well as the obtained reliability rate. Such value is defined
as the ratio between the number of optimizations in which the best indi-
vidual was obtained and the total number of optimizations.

%PENALTY.CONSTANT.FACTOR
<K>

(Value greater than 0)
Sets the penalty parameter when the Static penalty scheme is used. This
is usually a very large number.

%PENALTY.METHOD
'pen'

-
Static
Deb2000
Adaptive
-

Sets the desired penalty method. The Static method depends on a user-
specified penalty parameter, while the Deb2000 and Adaptive methods
are user-independent. Such method is used to evaluate the penalized
objective function of each trial project.

%POPULATION.SIZE
[PopSize]

(Value greater than 0)
Sets the size of the current population. Usually, increasing the
population size leads to better reliability rates. However, this can
lead to a high computational cost. When using the Island Model, each
island has PopSize individuals.

%PROBLEM.TYPE
'prob'

-
GSuite1
GSuite6
Colville8
Booth
Rastrigin
PVessel
Beam
-

Specifies the problem type which will be solved.

%SELECTION.METHOD
'sel'

-
Ranking
FitnessProportional
-

Specifies the selection method which will be used to determine
which individuals are going to participate in the crossover procedure.

-----------------------------------------------------------------------
READING THE OUTPUT FILE
-----------------------------------------------------------------------

The output file (.out extension) starts with a copy of the .opt input
data file, followed by the optimization results. Below is an explana-
tion of each result label:

%OPTIMIZATION NUMBER

Indicates to which optimization number the results pertain.

%RESULT.BEST.INDIVIDUALS

Prints the objective function for the best individual of the current
population for each generation.

%RESULT.AVERAGE.INDIVIDUALS

Prints the average objective function of the current population
for each generation.

%RESULT.INDIVIDUAL.EVALUATIONS

Indicates the number of times an individual was evaluated, which
also indicates how many times the analysis procedure was executed.

%RESULT.OBJECTIVE.FUNCTION

Prints the objective function of the best individual.

%RESULT.PENALIZED.OBJECTIVE.FUNCTION

Prints the penalized objective function of the best individual. For
a feasible individual, this has the same value as the objective
function.

%OPTIMIZATION.VARIABLES

Prints the design variables of the best individual in genotype
representation.

%PROBLEM.VARIABLES

Prints the design variables of the best individual in phenotype
representation, i.e. showing the actual value of each variable.

%CONSTRAINT.VALUES

Prints the value of each constraint in tabular form.
