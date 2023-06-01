def loadData(filename):
    with open(filename, 'r') as f:
        processor_count = int(f.readline().strip())
        process_count = int(f.readline().strip())
        processes = []
        for i in range(process_count):
            processes.append(int(f.readline().strip()))
    
    return {
        'processors': processor_count,
        'processes': processes
    }