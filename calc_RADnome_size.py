fin = open('working_data/test_run.RADnome.contig_positions.txt','rU')


total_bp = 0
for count, line in enumerate(fin):
    # if count > 10:
    #     break
    chrm, start, stop = line.strip().split(" ")
    start, stop = int(start), int(stop)
    size = stop - start
    if "NonOverlapping" in chrm:
        size = size - 50

    total_bp += size

print total_bp