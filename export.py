def export(filename, processors, execution_times, assignments):
    processor_usage = [0] * processors
    for i in range(len(assignments)):
        proc = assignments[i]
        processor_usage[proc] += execution_times[i]

    with open(f'export/{filename}.txt', 'w') as f:
        for usage in processor_usage:
            f.write(f'{usage}\n')