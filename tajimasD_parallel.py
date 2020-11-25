# function to calculate tajima's d from MSA
import Bio
from Bio import AlignIO #to parse the FASTA alignment
from Bio import SeqIO
import argparse
from multiprocessing import Pool, Lock, Process, Queue, current_process
import queue

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',dest='input',help='the FASTA multi-sequence alignment input')
parser.add_argument('-o','--output',dest='output',help='prefix of the csv file output')
parser.add_argument('-w','--window',dest='window',help='the window size to calculate tajima\'s D')
parser.add_argument('-s','--step',dest='step',help='the step size to instruct where to place the next window')
parser.add_argument('-t','--threads',dest='threads',help='number of cores to enable parallelization')
args = parser.parse_args()

#input fasta
align = AlignIO.read(args.input, "fasta")
out = open(args.output,"w+")

#functions for calculating tajimasD
def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

#do pairwise comparisons between the sequences
#collect them in a list and multiply by two
#divide the product by 'n choose 2'
def calculate_pi(alignment):
    align_len = len(align)
    counter = 0
    distances = []
    for i in range(align_len):
        j = 0
        j = j + counter
        while(j<align_len):
            if i == j:
                j = j + 1
                continue
            else:
                seq1 = alignment[i]
                seq2 = alignment[j]
                distances.append(hamming_distance(seq1,seq2))
                j = j + 1
        counter = counter + 1
    pi = (sum(distances)*2)/(align_len*(align_len-1))
    return pi

#calculate the number of segregating sites in the alignment
#then we calculate the value 'a' which is a summation of the
#divide the two value
def calculate_theta(alignment):
    seg_sites = 0
    for i in range(len(alignment[1].seq)):
        align_col = alignment[:,i]
        if len(list(set(align_col))) > 1:
            seg_sites = seg_sites + 1
        else:
            continue
    a = 0
    for i in range(1,len(alignment)):
        a = a + (1/i)
    return (seg_sites/a)

def tajimas_d(alignment):
    return (calculate_pi(alignment) - calculate_theta(alignment))

def tajima_run(task_to_accomplish):
    while True:
        try:
            '''
                try to get task from the queue. get_nowait() function will
                raise queue. Empty exception if the queue is empty.
                queue(False) function would do the same task also.
            '''
            region = task_to_accomplish.get_nowait()
        except queue.Empty:
            print("queue is empty")
            break
        else:
            '''
                if no exception has been raised, add the task completion
                message to task_that_are_done queue
            '''

            start = int(region.split(",")[0])
            end = int(region.split(",")[1])
            seq = align[:,start:end]
            tajimaD = tajimas_d(seq)
            out.write("{},{},{}".format(start,end,tajimaD))

    return True

def main():
    '''
    first we calculate how long the alignment is to break down the analysis into windows defined by the user.
    the number of processes are then evaluated based on what the user specifies.
    '''
    align_len = len(align[1])

    regions = []
    for i in range(0,align_len-(int(args.window)-1),int(args.step)):
        regions.append("{},{}".format(i,i+int(args.window)))


    number_of_task = len(regions)
    number_of_processes = int(args.threads)
    tasks_to_accomplish = Queue()
    tasks_that_are_done = Queue()
    processes = []

    for i in range(len(regions)):
        tasks_to_accomplish.put(regions[i])
    print(tasks_to_accomplish.qsize())

    # creating processes
    for w in range(number_of_processes):
        p = Process(target=tajima_run, args=(tasks_to_accomplish,))
        processes.append(p)
        p.start()

    # completing process
    for p in processes:
        p.join()

    return True

#run the main function where multiprocessing can occur
main()
out.close()
