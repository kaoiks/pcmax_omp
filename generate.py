from random import randint

with open('data.txt', 'w') as f:
    processor_count = randint(3, 25)
    process_count = randint(5 * processor_count, 250)

    processor_usage = [0] * processor_count
    processes = []

    for _ in range(process_count):
        duration = randint(1, 50)
        processes.append(duration)

        # Przyporządkuj do najmniej zajętego procesora
        proc_index_min = min(range(processor_count), key=processor_usage.__getitem__)
        processor_usage[proc_index_min] += duration

    max_usage = max(processor_usage)
    
    for i in range(processor_count):
        fill_length = max_usage - processor_usage[i]
        fill_length += randint(0, 1)
        processes.append(fill_length)
        processor_usage[i] += fill_length

    process_count = len(processes)
    cmax = max(processor_usage)

    print(f'Wygenerowano {process_count} procesów dla {processor_count} procesorów')
    print(f'Cmax* = {cmax}')

    f.write(f'{processor_count}\n{process_count}\n')
    for duration in processes:
        f.write(f'{duration}\n')