#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <climits>


struct Population {
    int cmax;
    std::vector<int> solution;
    
};

double GREEDY_MUTATION_CHANCE = 0.3;            // Prawdopodobieństwo, że mutacja będzie zachłanna
std::chrono::seconds MAX_DURATION = std::chrono::seconds(300);      // Maksymalny czas pracy w sekundach
double MAX_ITERATIONS = 1e6;                    //Maksymalna liczba iteracji
// double MAX_ITERATIONS = 50000;                    //Maksymalna liczba iteracji
double MIGRATION_CHANCE = 0.002;                //Prawdopodobieństwo migracji
double MUTATIONS_IN_SOLUTION = 2;               //Liczba mutacji w rozwiązaniu
int POPULATION_SIZE = 100;                      //Rozmiar populacji
double POPULATION_TO_DIE = 0.14;                //Odsetek populacji, który zginie w iteracji
double RANDOM_SOLUTIONS = 0.5;                  //Odsetek losowych rozwiązań w pierwotnej populacji
double SOLUTION_CROSSOVER_CHANCE = 0.85;        //Prawdopodobieństwo, że rozwiązanie się rozmnoży
double SOLUTION_MUTATION_CHANCE = 0.05;         //Prawdopodobieństwo, że w rozwiązaniu zajdzie mutacja

//Diagnostyka
bool ENABLE_EXPORT = false;
int PRINT_STATS_FREQ = 100;            //Co ile iteracji wyświetlać status

int all_time_best = INT_MAX;
int iterations;
Population best_solution;

#pragma omp threadprivate(iterations)

auto start_time = std::chrono::high_resolution_clock::now();

//Dane dla algorytmu
std::vector<int> execution_times = std::vector<int>();
int processor_count = 0;


//functions
std::vector<int> buildSolutionGreedy();
int measureSolutionCmax(std::vector<int> solution);
bool canContinue();
void performMigration(std::vector<Population>& population);
void doGeneticIteration(std::vector<Population>& population);
void performCrossOvers(std::vector<Population>& population);
void performMutations(std::vector<Population>& population);
void sortPopulation(std::vector<Population>& population);
std::vector<Population> crossOver(Population& parent1, Population& parent2);
void greedyMutate(Population& solution);
void mutate(Population& solution);


bool compareByCmax(const Population& p1, const Population& p2) {
    return p1.cmax < p2.cmax;
}

struct Data {
    int processors;
    std::vector<int> processes;
};

Data loadData(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    int process_count = 0;
    Data data;
    file >> data.processors;
    file >> process_count;
    
    for (int i = 0; i < process_count; i++) {
        int process;
        file >> process;
        // std::cout << process << std::endl; 
        data.processes.push_back(process);
    }
    
    file.close();

    return data;
}

//Generuje zestaw początkowych rozwiązań
std::vector<Population> generateInitialSolutions(int population_size){
    int process_count = execution_times.size();
    std::vector<Population>population(population_size);
    
    std::vector<int> solution(process_count, 0);
    for(int i = 0; i < population_size; ++i){
        double division_result = static_cast<double>(i)  / population_size;
        if (division_result < RANDOM_SOLUTIONS){
            for(int j = 0; j < processor_count; ++j){
                int range = processor_count - 1 - 0 + 1;
                int num = rand() % range + 0;
                solution[j] = num;
            }

        }
        else{
            solution = buildSolutionGreedy();
        }
        int cmax = measureSolutionCmax(solution);

        Population population_element;
        population_element.cmax = cmax;
        population_element.solution = solution;

        population[i] = population_element;
    }

    return population;
    
}


// Punkt wejściowy algorytmu
void genetic(int proc_count, std::vector<int> exec_times){
    processor_count = proc_count;
    execution_times = exec_times;

    start_time = std::chrono::high_resolution_clock::now();
    auto population = generateInitialSolutions(POPULATION_SIZE);
    std::sort(population.begin(), population.end(), compareByCmax);

    // best_cmaxes.push_back(population[0].cmax);
    // std::cout << population.size()<< std::endl;
    best_solution = population[0];
    all_time_best = population[0].cmax;
    while(iterations < MAX_ITERATIONS){
        doGeneticIteration(population);
        if(population[0].cmax < all_time_best){
            // if (population[0].cmax <= 9775){
            //     break;
            // }
            all_time_best = population[0].cmax;
            best_solution = population[0];
            std::cout << iterations << " " << all_time_best<< std::endl;
            // std::cout << iterations << " " << all_time_best<< std::endl;
        }
        // // std::cout << population[0].cmax<< std::endl;
        // best_cmaxes.push_back(population[0].cmax);
        // if (best_cmaxes.back() < all_time_best){
        //     all_time_best = best_cmaxes.back();
        //     best_solution = population[0];
        //     std::cout << iterations << " " << all_time_best<< std::endl;
        // }
    
        ++iterations;
       
    }
    std::cout << iterations << " " << all_time_best<< std::endl;
    auto count_duration = std::chrono::high_resolution_clock::now() - start_time;
    double duration_seconds = std::chrono::duration<double>(count_duration).count();
    std::cout << duration_seconds<< std::endl;
    // exit(0);
}


//Buduje rozwiązanie algorytmem zachłannym
std::vector<int> buildSolutionGreedy(){
    std::vector<int> solution(execution_times.size(), 0);
    std::vector<int> processor_usage(processor_count, 0);
    
    std::vector<int> order(execution_times.size(), 0);
    for(int i = 0; i < execution_times.size(); ++i){
        order[i] = i;
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(order.begin(), order.end(), gen);

    for(int i : order){
        auto min_element = std::min_element(processor_usage.begin(), processor_usage.end());
        if (min_element != processor_usage.end()) {
            int proc_index_min = std::distance(processor_usage.begin(), min_element);
            processor_usage[proc_index_min] += execution_times[i];
            solution[i] = proc_index_min;
        }

    }
    return solution;
}


void doGeneticIteration(std::vector<Population>& population){
    performMigration(population);
    performCrossOvers(population);  
    performMutations(population);
    sortPopulation(population);
}


void sortPopulation(std::vector<Population>& population){
    std::sort(population.begin(), population.end(), [](const Population& s1, const Population& s2) {
        return s1.cmax < s2.cmax;
    });

}

void performMigration(std::vector<Population>& population){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    if (dist(gen) <= MIGRATION_CHANCE){
        if (all_time_best < population[0].cmax){
            population[population.size()-1] = best_solution;
        }
        else{
            std::vector<int> solution = buildSolutionGreedy();
            Population new_population;
            new_population.cmax = measureSolutionCmax(solution);
            new_population.solution = solution;
            population[population.size()-1] = new_population;
        }
    } 
}

void performCrossOvers(std::vector<Population>& population){
    float cross_overs = population.size() * POPULATION_TO_DIE;
    // std::vector<Population> (population.size());
    //zamiast 15 losowych umiera 15 najgorszych py vs cpp
    // Wypełnij tablicę dziećmi
    int idx = population.size() - cross_overs;
    // while (idx < cross_overs && idx + 1 < population.size()){
    while (idx < population.size()){

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(0, population.size() - 1);

        std::uniform_real_distribution<> dist2(0.0, 1.0);
        if (dist2(gen) >= SOLUTION_CROSSOVER_CHANCE){
            continue;
        }

        int parent1 = dist(gen);
        int parent2 = dist(gen);
        while(parent1 == parent2){
            parent2 = dist(gen);
        }
        
        
        std::vector<Population> child_arr = crossOver(population[parent1], population[parent2]);
        for (auto child : child_arr){
            population[idx] = child;
            idx += 1;
        }
    }
    
}


std::vector<Population> crossOver(Population& parent1, Population& parent2){
    std::vector<int>& solution1 = parent1.solution;
    std::vector<int>& solution2 = parent2.solution;

    std::vector<int> offspring_1(solution1);
    std::vector<int> offspring_2(solution2);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, solution1.size() - 1);
    int start = dist(gen);
    std::uniform_int_distribution<> dist2(start, solution1.size());
    int end = dist2(gen);

    for (int i = start; i < end; ++i) {
        offspring_1[i] = solution2[i];
        offspring_2[i] = solution1[i];
    }

    Population child1;
    child1.cmax = measureSolutionCmax(offspring_1);
    child1.solution = offspring_1;
    

    Population child2;
    child2.cmax = measureSolutionCmax(offspring_2);
    child2.solution = offspring_2;
    std::vector<Population> children;
    children.push_back(child1);
    children.push_back(child2);

    return children;
}


void performMutations(std::vector<Population>& population){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    
    for(int i = 0; i < population.size(); ++i){
        if(dist(gen) >= SOLUTION_MUTATION_CHANCE){
            continue;
        }
        if (dist(gen) <= GREEDY_MUTATION_CHANCE){
            greedyMutate(population[i]);
        }
        else{
            mutate(population[i]);
        }
    }
}

void mutate(Population& solution){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, solution.solution.size() - 1);

    for(int i = 0; i < MUTATIONS_IN_SOLUTION; ++i){
        int pos1 = dist(gen);
        int pos2 = dist(gen);
        int solution_pos1 = solution.solution[pos1];
        int solution_pos2 = solution.solution[pos2];

        solution.solution[pos1] = solution_pos2;
        solution.solution[pos2] = solution_pos1;

    }
    solution.cmax = measureSolutionCmax(solution.solution);
}


void greedyMutate(Population& solution){
    std::vector<int> processor_usage(processor_count, 0);
    // std::vector<int> new_assignment = solution.solution;

    for(int i = 0; i < solution.solution.size(); ++i){
        int proc = solution.solution[i];
        processor_usage[proc] += execution_times[i];
    }

    for(int j = 0; j < 5; ++j){
        int min_load[2] = {INT_MAX, -1};
        int max_load[2] = {0, -1};
        for(int i = 0; i < processor_usage.size(); ++i){
            if (processor_usage[i] < min_load[0]){
                min_load[0] = processor_usage[i];
                min_load[1] = i;
            }
            if (processor_usage[i] > max_load[0]){
                max_load[0] = processor_usage[i];
                max_load[1] = i;
            }
        }
        
        int shortest_task[2] = {INT_MAX, -1};
        for(int i = 0; i < solution.solution.size(); ++i){
            if (solution.solution[i] == max_load[1] && execution_times[i] < shortest_task[0]){
                shortest_task[0] = execution_times[i];
                shortest_task[1] = i;
            }
        }

        solution.solution[shortest_task[1]] = min_load[1];
        processor_usage[min_load[1]] += shortest_task[0];
        processor_usage[max_load[1]] -= shortest_task[0];
    }
    // Population new_population;
    // new_population.cmax = measureSolutionCmax(new_assignment);
    // new_population.solution = new_assignment;
    // return new_population;
    solution.cmax =  measureSolutionCmax(solution.solution);
}


// Mierzy Cmax rozwiązania (im mniej tym lepiej)
int measureSolutionCmax(std::vector<int> solution){
    std::vector<int> processor_occupancy(processor_count, 0);
    for(int i = 0; i < solution.size(); ++i){
        int proc = solution[i];
        processor_occupancy[proc] += execution_times[i];
    }
    auto max_element = std::max_element(processor_occupancy.begin(), processor_occupancy.end());
    int max_value = static_cast<int>(*max_element);

    return max_value;
}


int main(int argc, char* argv[]){
    srand(time(0));
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
        return 1; 
    }

    std::string fname = argv[1]; 

    Data data = loadData(fname); 
    #pragma omp parallel num_threads(8) shared(all_time_best, best_solution)
    {
        Data data = loadData(fname);
        iterations = 0;
        genetic(data.processors, data.processes);
    }
    return 0;
}
   