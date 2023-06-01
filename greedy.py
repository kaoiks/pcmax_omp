from load import loadData
import sys

def greedy(processors: int, processes: list):
    tab = []
    sums = []
    
    #tablica procesorow
    for i in range(0,processors):
        tab.append([])
        sums.append(0)
    
    for process in processes:
        min=sums[0]
        minIndex=0
        for i in range(1,processors):
            if sums[i] < min:
                min=sums[i]
                minIndex = i
        tab[minIndex].append(process)
        sums[minIndex]+=process
    return tab


def main():
    fname = 'data.txt'
    if len(sys.argv) >= 2:
        fname = sys.argv[1]

    data = loadData(fname)

    result = greedy(data['processors'], data['processes'])
    for i in range(len(result)):
        print(f'Procesor {i+1}: '+str(result[i]))
    
    sums = []
    for res in result:
        sums.append(sum(res))
    print('Cmax: ' + str(max(sums)))

main()