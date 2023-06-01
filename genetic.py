from export import export
from math import inf, log10
from load import loadData
from random import random, randint, shuffle
import sys
from time import time, sleep
from typing import List, Tuple

#
# Populacja przechowywana jest jako lista par. Pierwszy element pary
# to Cmax rozwiązania, a drugi - tablica przypisująca zadanie do procesora.
# Ta realizacja pozwala na oszczędność w czasie działania i jest o ok. 70% szybsza
# niż przechowywanie poszczególnych rozwiązań jako obiektów klasy Solution.
#

# Konfiguracja algorytmu
GREEDY_MUTATION_CHANCE = 0.3            # Prawdopodobieństwo, że mutacja będzie zachłanna
MAX_DURATION = 300                      # Maksymalny czas pracy w sekundach
MAX_ITERATIONS = 1e6                    # Maksymalna liczba iteracji
MIGRATION_CHANCE = 0.002                # Prawdopodobieństwo migracji
MUTATIONS_IN_SOLUTION = 2               # Liczba mutacji w rozwiązaniu
POPULATION_SIZE = 100                   # Rozmiar populacji
POPULATION_TO_DIE = 0.15                # Odsetek populacji, który zginie w iteracji
RANDOM_SOLUTIONS = 0.5                  # Odsetek losowych rozwiązań w pierwotnej populacji
SOLUTION_CROSSOVER_CHANCE = 0.85        # Prawdopodobieństwo, że rozwiązanie się rozmnoży
SOLUTION_MUTATION_CHANCE = 0.05         # Prawdopodobieństwo, że w rozwiązaniu zajdzie mutacja

# Diagnostyka
ENABLE_EXPORT = False
PRINT_STATS_FREQ = 100                  # Co ile iteracji wyświetlać status

# Zmienne związane z pracą programu
all_time_best = inf
best_cmaxes = []
best_solution = None
iterations = 0
start_time = 0.0

# Dane dla algorytmu
execution_times = []
processor_count = 0

# Punkt wejściowy algorytmu
def genetic(proc_count: int, exec_times: List[int]) -> None:
    global all_time_best, best_cmaxes, best_solution, iterations, start_time
    global execution_times, processor_count

    execution_times = exec_times
    processor_count = proc_count
    start_time = time()

    population = generateInitialSolutions(POPULATION_SIZE)
    sortPopulation(population)

    best_cmaxes.append(population[0][0])
    best_solution = population[0]
    while canContinue():
        if iterations % PRINT_STATS_FREQ == 0 and ENABLE_EXPORT:
            export(iterations, processor_count, execution_times, best_solution[1])
        population = doGeneticIteration(population)
        best_cmaxes.append(population[0][0])
        if best_cmaxes[-1] < all_time_best:
            all_time_best = best_cmaxes[-1]
            best_solution = population[0]
        iterations += 1
        printStats()
    if ENABLE_EXPORT:
        export(iterations, processor_count, execution_times, best_solution[1])


# Sprawdza czas trwania, Cmax rozwiązania i ew. inne metryki i decyduje czy kontynuować
def canContinue() -> bool:
    global best_cmaxes, iterations, start_time
    global MAX_DURATION, MAX_ITERATIONS
    duration = time() - start_time

    # Po zakończeniu wykonywania programu można podnieść limity
    if iterations >= MAX_ITERATIONS:
        print('[?] Iteration limit reached. Continue for next 50? [y/N]? ', end='')
        answer = input()
        if answer.lower() == 'y':
            MAX_ITERATIONS += 50
    
    if duration >= MAX_DURATION:
        print('[?] Time limit reached. Continue for next 60 seconds? [y/N]? ', end='')
        answer = input()
        if answer.lower() == 'y':
            MAX_DURATION = (time() + 60) - start_time

    return (iterations < MAX_ITERATIONS
        and duration <= MAX_DURATION)


# Wykonuje iterację algorytmu genetycznego
def doGeneticIteration(population: List[Tuple[int, List[int]]]) -> List[Tuple[int, List[int]]]:
    population = performCrossOvers(population)
    performMutations(population)
    performMigration(population)
    sortPopulation(population)
    return population


# Generuje zestaw początkowych rozwiązań
def generateInitialSolutions(population_size: int) -> List[Tuple[int, List[int]]]:
    global execution_times, processor_count
    process_count = len(execution_times)
    population = [None] * population_size
    for i in range(population_size):
        if i / population_size < RANDOM_SOLUTIONS:
            solution = [ randint(0, processor_count-1) for _ in range(process_count) ]
        else:
            solution = buildSolutionGreedy()
        cmax = measureSolutionCmax(solution)
        population[i] = (cmax, solution)
    return population


# Buduje rozwiązanie algorytmem zachłannym
def buildSolutionGreedy() -> List[int]:
    global execution_times, processor_count
    solution = [0] * len(execution_times)
    processor_usage = [0] * processor_count
    order = list(range(len(execution_times)))
    
    shuffle(order)
    
    for i in order:
        #print(range(len(processor_usage)))
        proc_index_min = min(range(len(processor_usage)), key=processor_usage.__getitem__)
        processor_usage[proc_index_min] += execution_times[i]
        solution[i] = proc_index_min
    return solution


# Wybiera rozwiązania do rozmnożenia i tworzy nową populację
# Jeśli dzieci jest mniej niż rodziców, to nowa populacja jest
# uzupełniana starą
def performCrossOvers(population: List[Tuple[int, List[int]]]) -> List[Tuple[int, List[int]]]:
    cross_overs = len(population) * POPULATION_TO_DIE
    children = [None] * len(population)

    # Wypełnij tablicę dziećmi
    idx = 0
    while idx < cross_overs and idx+1 < len(children):
        # Wybierz dwa różne rozwiązania do rozmnożenia
        parent1 = select(population)
        parent2 = select(population)
        while parent1 == parent2:
            parent2 = select(population)

        if random() >= SOLUTION_CROSSOVER_CHANCE:
            continue
        
        child_arr = crossOver(population[parent1], population[parent2])
        for child in child_arr:
            children[idx] = child
            idx += 1

    # Uzupełnij populację osobnikami z poprzedniej
    parents_to_survive = set()
    while idx < len(children):
        parent = select(population)
        if parent in parents_to_survive:
            continue
        parents_to_survive.add(parent)
        children[idx] = population[parent]
        idx += 1

    return children


# Rozmnaża rozwiązania
def crossOver(parent1: Tuple[int, List[int]], parent2: Tuple[int, List[int]]) -> List[Tuple[int, List[int]]]:
    solution1 = parent1[1]
    solution2 = parent2[1]
    offspring_1 = solution1.copy()
    offspring_2 = solution2.copy()

    start = randint(0, len(solution1) - 1)
    end = randint(start, len(solution1))

    offspring_1[start:end] = solution2[start:end]
    offspring_2[start:end] = solution1[start:end]
    # print(offspring_1)
    
    child1 = (measureSolutionCmax(offspring_1), offspring_1)
    child2 = (measureSolutionCmax(offspring_2), offspring_2)
    return [child1, child2]


# Wybiera rozwiązania i dokonuje mutacji
def performMutations(population: List[Tuple[int, List[int]]]) -> None:
    global execution_times, processor_count
    for i in range(len(population)):
        if random() >= SOLUTION_MUTATION_CHANCE:
            continue
        if random() <= GREEDY_MUTATION_CHANCE:
            population[i] = greedyMutate(population[i])
        else:
            population[i] = mutate(population[i])


# Mutuje rozwiązanie i zwraca nową kopię
def mutate(solution: Tuple[int, List[int]]) -> Tuple[int, List[int]]:
    global execution_times, processor_count
    new_solution = solution[1].copy()
    for _ in range(MUTATIONS_IN_SOLUTION):
        pos1 = randint(0, len(new_solution) - 1)
        pos2 = randint(0, len(new_solution) - 1)
        new_solution[pos1] = solution[1][pos2]
        new_solution[pos2] = solution[1][pos1]

    return (measureSolutionCmax(new_solution), new_solution)


def greedyMutate(solution: Tuple[int, List[int]]) -> Tuple[int, List[int]]:
    global execution_times, processor_count
    processor_usage = [0] * processor_count

    new_assignment = solution[1].copy()

    for i in range(len(solution[1])):
        proc = solution[1][i]
        processor_usage[proc] += execution_times[i]

    for _ in range(5):
        min_load = [inf, -1]    # [value, index]
        max_load = [0, -1]
        for i in range(len(processor_usage)):
            if processor_usage[i] < min_load[0]:
                min_load[0] = processor_usage[i]
                min_load[1] = i
            if processor_usage[i] > max_load[0]:
                max_load[0] = processor_usage[i]
                max_load[1] = i
        shortest_task = [inf, -1]    # [value, index]
        for i in range(len(new_assignment)):
            if new_assignment[i] == max_load[1] and execution_times[i] < shortest_task[0]:
                shortest_task[0] = execution_times[i]
                shortest_task[1] = i
        new_assignment[shortest_task[1]] = min_load[1]
        processor_usage[min_load[1]] += shortest_task[0]
        processor_usage[max_load[1]] -= shortest_task[0]

    return (measureSolutionCmax(new_assignment), new_assignment)


# Sortuje populację od najlepszych rozwiązań
def sortPopulation(population: List[Tuple[int, List[int]]]) -> None:
    population.sort(key=lambda s: s[0])


# Dokonuje migracji dobrego rozwiązania do populacji
def performMigration(population: List[Tuple[int, List[int]]]) -> None:
    if random() <= MIGRATION_CHANCE:
        if all_time_best < best_cmaxes[-1]:
            population[-1] = best_solution
        else:
            solution = buildSolutionGreedy()
            population[-1] = (measureSolutionCmax(solution), solution)


# Mierzy Cmax rozwiązania (im mniej tym lepiej)
def measureSolutionCmax(solution: List[int]) -> None:
    global execution_times, processor_count
    processor_occupancy = [0] * processor_count
    for i in range(len(solution)):
        proc = solution[i]
        processor_occupancy[proc] += execution_times[i]
    return max(processor_occupancy)


# Zwraca indeks rozwiązania z puli
def select(population: List[Tuple[int, List[int]]]) -> int:
    return randint(0, len(population) - 1)


# Wypisuje statystyki co określoną liczbę iteracji
def printStats() -> None:
    if iterations % PRINT_STATS_FREQ != 0:
        return
    duration = round(time() - start_time, 1)
    digits = int(log10(MAX_ITERATIONS)) + 1
    print(f'[{iterations: >{digits}}]: {duration: >6}s elapsed, Cmax: {best_cmaxes[-1]}, ATB: {all_time_best}')


# Wypisuje końcowe statystyki
def printFinalStats() -> None:
    global all_time_best, best_cmaxes, iterations, start_time
    global execution_times, processor_count
    optimum = round(sum(execution_times) / processor_count, 1)

    duration = time() - start_time
    duration = round(duration, 1)
    print('\nFinished job')
    print(f'    {iterations} iterations performed')
    print(f'    {duration} seconds')
    print(f'    \033[1;32m{all_time_best} = Cmax\033[0m')
    print(f'    {optimum} = Cmax* <divisible tasks>')
    print(f'    {best_cmaxes[0]} -> {all_time_best} genetic progress')


def main() -> None:
    fname = 'm50n200.txt'
    if len(sys.argv) >= 2:
        fname = sys.argv[1]

    data = loadData(fname)
    genetic(data['processors'], data['processes'])
    printFinalStats()

main()