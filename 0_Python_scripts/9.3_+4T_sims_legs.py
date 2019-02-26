## IMPORTS ###

import os
import csv
import numpy as np
import time
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '2_LEGs_fasta'

### FUNCTIONS ### - Get main to work without parallelizing, then edit to include the hashed out commands

def main():
    filenames = get_files(source)
    workers = 10
    processes = run_in_parallel(filenames, ['foo', source], folder_parse, workers=workers)

    csv_total = []

    for process in processes:
        output = process.get()
        csv_total.extend(output)

    write_outputs(csv_total)


def folder_parse(filenames, source):

    csv_total = []

    for i,f in enumerate(filenames):

        print ('Doing genome {0}/{1}'.format(i+1, len(filenames)))

        path = os.path.join(source, f)
        raw = open(path).read()

        if '>' in raw:

            #Split into gene list for use in functions
            genes = raw.strip().split('>')
            genes_c = list(filter(None, genes))

            #Obtain accession
            a = raw.split('>')
            a1 = a[1].split(';')
            accession = a1[0]

            #Obtain GC and GC3
            a2 = a1[1].strip()
            a3 = a2.split("=")
            GC = float(a3[1])

            a4 = a1[2].strip()
            a5 = a4.split("=")
            GC3 = float(a5[1])

            #Get a list of all genes containing +4T, then calculate additional stop codon frequency at +1 to +6
            fourth_t = get_fourth_t(genes_c)

            #Only consider genomes where there are 4th T genes
            if len(fourth_t) > 0:

                n, real1, real2, real3, real4, real5, real6 = get_counts(fourth_t)
                rfreq1, rfreq2, rfreq3, rfreq4, rfreq5, rfreq6 = get_freqs(n, real1, real2, real3, real4, real5, real6)

                #Combine UTR sequences for later use in calculating dinucleotide content
                total_utr = get_total_UTR(genes_c)

                #Simulate 10,000 null sequences, get frequencies
                sim_list = get_simulations(total_utr)
                n_sims, exp1, exp2, exp3, exp4, exp5, exp6 = get_sim_counts(sim_list)
                sfreq1, sfreq2, sfreq3, sfreq4, sfreq5, sfreq6 = get_freqs(n_sims, exp1, exp2, exp3, exp4, exp5, exp6)

                #Create lists for stats
                raw_ob_counts = [real1, real2, real3, real4, real5, real6]
                total_utrs = [n, n, n, n, n, n]
                ex_frequencies = [sfreq1, sfreq2, sfreq3, sfreq4, sfreq5, sfreq6]

                #Create nested list to iterate through
                stats_list = [raw_ob_counts, total_utrs, ex_frequencies]
                nested = [list(i) for i in zip(*stats_list)]

                #Calculate and add observed frequency, standard deviation, z-score and binomial p-value to the nested list
                for i in nested:
                    observed_freq = i[0] / i[1]
                    i.append(observed_freq)

                    standard_dev = (10000 * (i[2]) * (1-i[2])) ** 0.5
                    i.append(standard_dev)

                    z_score = (i[0] - (i[1] * i[2])) / i[4]
                    i.append(z_score)

                    p = st.binom_test(i[0], i[1], i[2])
                    i.append(p)

                #Flatten list
                output_flattened = [item for sublist in nested for item in sublist]

                #Add Accession and GC content information to the list, which will be added to the output CSV
                csv_nested = [[accession], output_flattened, [GC], [GC3]]
                csv_list = [item for sublist in csv_nested for item in sublist]
                csv_total.append(csv_list)

    return csv_total


def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example, a list of gene names)
    args: a list of arguments to the function. Put in "foo" in place of the argument you are parallelizing over.
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for elem in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = elem
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()

    return(results)


def get_files(source):

    ''' Get list of file names to iterate through '''

    files = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            files.append(f)

    return files


def get_fourth_t(chunks):

    fourth_t = []

    for i in chunks:

        #Obtain UTRs
        split_genes = i.split('\n')
        utr = split_genes[1]

        if utr[3] == 't':
            fourth_t.append(utr)

    return fourth_t


def get_counts(utr_list):

    #Start counters for each position
    pos1 = 0
    pos2 = 0
    pos3 = 0
    pos4 = 0
    pos5 = 0
    pos6 = 0

    #Calculate length
    n = len(utr_list)

    #Add to counters depending on additional stop codon presence/absence
    for utr in utr_list:

        #Convert into codons
        codon_seq = [utr[i:i+3] for i in range(0, len(utr), 3)]

        if codon_seq[1] == 'taa' or codon_seq[1] == 'tag' or codon_seq[1] == 'tga':
            pos1 += 1
        if codon_seq[2] == 'taa' or codon_seq[2] == 'tag' or codon_seq[2] == 'tga':
            pos2 += 1
        if codon_seq[3] == 'taa' or codon_seq[3] == 'tag' or codon_seq[3] == 'tga':
            pos3 += 1
        if codon_seq[4] == 'taa' or codon_seq[4] == 'tag' or codon_seq[4] == 'tga':
            pos4 += 1
        if codon_seq[5] == 'taa' or codon_seq[5] == 'tag' or codon_seq[5] == 'tga':
            pos5 += 1
        if codon_seq[6] == 'taa' or codon_seq[6] == 'tag' or codon_seq[6] == 'tga':
            pos6 += 1

    return n, pos1, pos2, pos3, pos4, pos5, pos6


def get_sim_counts(utr_list):

    #Start counters for each position
    pos1 = 0
    pos2 = 0
    pos3 = 0
    pos4 = 0
    pos5 = 0
    pos6 = 0

    #Calculate length
    n = len(utr_list)

    #Add to counters depending on additional stop codon presence/absence
    for utr in utr_list:

        #Convert into codons
        codon_seq = [utr[i:i+3] for i in range(0, len(utr), 3)]

        if codon_seq[0] == 'taa' or codon_seq[0] == 'tag' or codon_seq[0] == 'tga':
            pos1 += 1
        if codon_seq[1] == 'taa' or codon_seq[1] == 'tag' or codon_seq[1] == 'tga':
            pos2 += 1
        if codon_seq[2] == 'taa' or codon_seq[2] == 'tag' or codon_seq[2] == 'tga':
            pos3 += 1
        if codon_seq[3] == 'taa' or codon_seq[3] == 'tag' or codon_seq[3] == 'tga':
            pos4 += 1
        if codon_seq[4] == 'taa' or codon_seq[4] == 'tag' or codon_seq[4] == 'tga':
            pos5 += 1
        if codon_seq[5] == 'taa' or codon_seq[5] == 'tag' or codon_seq[5] == 'tga':
            pos6 += 1

    return n, pos1, pos2, pos3, pos4, pos5, pos6


def get_freqs(n, count1, count2, count3, count4, count5, count6):

    freq1 = count1 / n
    freq2 = count2 / n
    freq3 = count3 / n
    freq4 = count4 / n
    freq5 = count5 / n
    freq6 = count6 / n

    return freq1, freq2, freq3, freq4, freq5, freq6


def get_total_UTR(chunks):

    total_utr = ''

    for i in chunks:

        #Obtain UTRs and add to total utr string
        split_genes = i.split('\n')
        utr = split_genes[1]
        total_utr += utr

    return total_utr


def get_simulations(total_utr):

    #No need to calculate probability of first base, as this will be selected to be T

    #Probability of next base - create dictionary with pool of next nuc options
    trans = {}
    for i in range(len(total_utr)-1):

        dinuc = total_utr[i:i+2]
        first_dinuc = dinuc[0]
        sec_dinuc = dinuc[1]

        if first_dinuc not in trans:
            trans[first_dinuc] = [sec_dinuc]
        else:
            trans[first_dinuc] += sec_dinuc

    #Build simulations
    total_list = []
    count = 0

    #This while loop determines how many simulations will be completed
    while count < 10000:

        codon_list = []
        sim = []

        np.random.seed()

        #First base must be T
        sim.append('t')

        #For one gene simulation
        while (len(sim) < 21):

            #Generate next base
            prev_base = sim[-1]
            next_base = np.random.choice(trans[prev_base], replace=True)
            sim.append(next_base)

        sim = "".join(sim)
        total_list.append(sim)
        count += 1

    return total_list


def write_outputs(csv_total):

    #Set headers for the CSV
    headers = ["Accession", "P1_observed", "P1_Total_UTRs", "P1_Expected_freq", "P1_Observed_freq", "P1_Standard_Deviation", "P1_Zscore", "P1__BinomPval",
        "P2_observed", "P2_Total_UTRs", "P2_Expected_freq", "P2_Observed_freq", "P2_Standard_Deviation", "P2_Zscore", "P2__BinomPval",
        "P3_observed", "P3_Total_UTRs", "P3_Expected_freq", "P3_Observed_freq", "P3_Standard_Deviation", "P3_Zscore", "P3__BinomPval",
        "P4_observed", "P4_Total_UTRs", "P4_Expected_freq", "P4_Observed_freq", "P4_Standard_Deviation", "P4_Zscore", "P4__BinomPval",
        "P5_observed", "P5_Total_UTRs", "P5_Expected_freq", "P5_Observed_freq", "P5_Standard_Deviation", "P5_Zscore", "P5__BinomPval",
        "P6_observed", "P6_Total_UTRs", "P6_Expected_freq", "P6_Observed_freq", "P6_Standard_Deviation", "P6_Zscore", "P6__BinomPval",
        "Genomic_GC", "Genomic_GC3"]

    create_csv(headers, csv_total)


def create_csv(headers, csv_total):

    filename = "9.1_+4T_sims_LEGs.csv"
    subdir = "4_Outputs/CSVs"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
