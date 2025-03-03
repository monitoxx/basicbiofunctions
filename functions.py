import re
import random
import pandas as pd


def gc_content(sequence):
    '''It processes a nucleotide sequence counting only Guanine and Cytosine bases
     for calculating the GC percentage of the sequence and additionally the total
     length. Returns GC content and Length in this particular order.  Refer to the example: www.github.com/monitoxx/Retrovirus-Analyzer-V2'''

    sequence = sequence.rstrip("\n")
    sequence = sequence.upper()
    G_count = sequence.count("G")
    C_count = sequence.count("C")
    GC_content = G_count + C_count
    GC_content = GC_content / len(sequence)
    GC_content = round(100 * GC_content, 2)
    return float(GC_content), int(len(sequence))
    pass

#print(gc_content.__doc__)


def sequence_reader(file_name, gene_name):
    '''Identifies a gene of interest and extracts the correspondent sequence.
    The input has to be the file name with extension, as it is on the file directory (Example: sequence.txt
    or sequence.fasta) and the gene name exactly as in the sequence (Example: [gene=gag]).
    The files have to be in coding sequences format and the gene in str. The output from this function is
    the Genebank code of the sequence and the sequence respectively. Refer to the example: www.github.com/monitoxx/Retrovirus-Analyzer-V2'''
    sequence = ''
    sequence_name = None

    try:
        with open(file_name, 'r') as file:
            include_sequence = False  # Flag for including the sequence
            for line in file:
                #include_sequence = False # Flag for including the sequence
                if line.startswith('>'):
                    # If header contains our gene of interest, activates the flag
                    include_sequence = gene_name in line
                    if include_sequence:
                        # Searches for the GenBank code
                        match = re.search(r'\|(.+?)\_', line)
                        if match:
                            sequence_name = match.group(1)  # Extracts code of Genbank
                elif include_sequence:  #Acumulates lines of sequence only when our flag is activated
                    line_ = line.strip()
                    sequence += line.strip()

    except FileNotFoundError:
        print(f"File {file_name} not found.")


    return sequence_name, sequence

def gc_list_generator(sequence, sequence_name):
    '''Generates a Guanine-Cytosine percentage list from sections of the nucleotide sequence entered. It does the same
    for a random shuffled sequence generated from the latter. The step length is 70 positions. It returns in the following
    order: the position list (for plotting ease purposes), the GC list, the random GC list and the GenBank code list
    from the sequence entered. Refer to the example: www.github.com/monitoxx/Retrovirus-Analyzer-V2'''
    position_list = []
    gc_list = []
    gc_list_random = []
    code_list = []

    # Generates a random sequences once
    seq_random_list = list(sequence)
    random.shuffle(seq_random_list)
    seq_random = ''.join(seq_random_list)

    for x in range(0, len(sequence), 70):
        # Fragments of the original and random sequences
        temp = sequence[x:x + 70]
        temp2 = seq_random[x:x + 70]

        # Calculates GC content
        gc_temp, _ = gc_content(temp)
        gc_temp2, _ = gc_content(temp2)

        # Adds dats to lists
        position_list.append(x)
        gc_list.append(gc_temp)
        gc_list_random.append(gc_temp2)
        code_list.append(sequence_name)

    return position_list, gc_list, gc_list_random, code_list

def df_creator(file_names, gene_name):
    '''A powerful function made for easier analysis, interpretation and plotting, of nucleotide sequences that want to be
    analyzed in a specific gene, in order to show relatedness between them or not. Initially made for comparison of different HIV
    sequences in the 'gag' gene.
     The input is a list of the file names of the sequences in CDS format in the
     directory (Example: file_names = ['sequence.txt', 'sequence2.txt', 'sequence3.txt', 'sequence4.txt']), and the gene name
     exactly like it is found in the CDS file of the sequence that wants to be analyzed (Example: '[gene=gag]').
    There are 6 outputs expected, they are, in order:
      df_sequences: a dataframe organized with the genbank code,
      the complete sequence in 1 line and the length of each one of them

      length: the minimum length sequence of them all, useful for plotting them in the same dimension and operating 'for'
      cicles with the correct index.

      df_gc_content: dataframe derived from df_sequences, indicates the genbank code, %GC content and length of the sequences

      df_gc_content_lists: dataframe with genbank codes, %GC content for every 70 nucleotides (the size of the step) and
      the positions for easy plotting of the results, and also the GC content for a randomized version of the nucleotide sequence.

      df_gc_content_lists_exploded: dataframe derived from df_gc_content_lists but expanded as columns and with additions
      such as quartiles for doing distributions.

      heatmap_data: a dataframe specifically for doing a heatmap with the data of the GC content from the original
      sequence, partitioned every 70 nucleotides. The pivoted version of df_gc_content_lists_exploded only with genbank code,
      position and GC content.
     '''
    # Create an empty dataFrame
    df_sequences = pd.DataFrame(columns=['GenBank Code', 'Sequence', 'Length'])
    df_gc_content = pd.DataFrame(columns=['GenBank Code', '% GC Content', 'Length'])
    df_gc_content_lists = pd.DataFrame(columns=['GenBank Code', '% GC Content', '% GC Content Random', 'Position'])

    # Iterate the files and add the data to the dataFrames
    for file_name in file_names:
        sequence_name, sequence = sequence_reader(file_name, gene_name)
        if sequence_name and sequence:  # Only adds up data if they are valid
            var_gc_content, len_sequence = gc_content(sequence)
            position_list, gc_list, gc_list_random, code_list = gc_list_generator(sequence, sequence_name)

            df_sequences = pd.concat([df_sequences, pd.DataFrame(
                {'GenBank Code': [sequence_name],
                 'Sequence': [sequence],
                 'Length': [len_sequence]
                 })], ignore_index=True)
            df_gc_content = pd.concat([df_gc_content if not df_gc_content.empty else None, pd.DataFrame(
                {'GenBank Code': [sequence_name],
                 '% GC Content': [var_gc_content],
                 'Length': [len_sequence]
                 })], ignore_index=True)
            df_gc_content_lists = pd.concat([df_gc_content_lists, pd.DataFrame(
                {'GenBank Code': [code_list],
                 '% GC Content': [gc_list],
                 '% GC Content Random': [gc_list_random],
                 'Position': [position_list]
                 })], ignore_index=True)

    #This will later help me for my spearman analysis, so that indexes match
    length = min(df_sequences['Length'])


    # Unfolds lists in lines
    df_gc_content_lists_exploded = df_gc_content_lists.explode(['% GC Content', '% GC Content Random', 'Position', 'GenBank Code'])
    # Makes sure that numeric columns are in the correct format:
    df_gc_content_lists_exploded['% GC Content'] = pd.to_numeric(df_gc_content_lists_exploded['% GC Content'])
    df_gc_content_lists_exploded['% GC Content Random'] = pd.to_numeric(df_gc_content_lists_exploded['% GC Content Random'])
    df_gc_content_lists_exploded['Position'] = pd.to_numeric(df_gc_content_lists_exploded['Position'])


    # Group data by GenBank Code
    grouped_data = df_gc_content_lists_exploded.groupby('GenBank Code')

    # Calculate quartiles for each group
    df_gc_content_lists_exploded['Quartile'] = grouped_data['% GC Content'].transform(pd.qcut, q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    # Repeat for random GC content
    df_gc_content_lists_exploded['Quartile Random'] = grouped_data['% GC Content Random'].transform(pd.qcut, q=4, labels=["Q1", "Q2", "Q3", "Q4"])

    # Creates a pivot matrix for heatmap:
    heatmap_data = df_gc_content_lists_exploded.pivot_table(
        index='GenBank Code', columns='Position', values='% GC Content', aggfunc='mean')


    #print(df_gc_content_lists_exploded)

    return df_sequences, length, df_gc_content, df_gc_content_lists, df_gc_content_lists_exploded, heatmap_data



def database_sample(file_name, gene_name):
    gene_name = str(gene_name)
    file_names = []
    try:
        with open(file_name, 'r') as file:
            sequence = ''
            sequence_name = None
            count = 0  # Contador de códigos detectados
            file_index = 1  # Índice de archivo

            while count <= 9:  # Se detiene después de 5 códigos detectados
                line = file.readline()
                if not line:
                    break  # Detener si se acaba el archivo

                if line.startswith('>'):
                    if sequence_name and sequence:
                        # Guardar la secuencia anterior en un archivo
                        file_name = f'output_{file_index}.txt'
                        file_names.append(file_name)
                        with open(file_name, 'w') as out_file:
                            out_file.write(f">|{sequence_name}_ [gene={gene_name}]\n{sequence}\n")

                        file_index += 1  # Incrementar el índice de archivo
                        sequence = ''  # Reiniciar secuencia

                    match = re.search(r'>(.+?)\:', line)
                    if match:
                        sequence_name = match.group(1)  # Extrae código de GenBank
                        count += 1
                else:
                    sequence += line.strip()  # Acumula secuencia

            # Guardar la última secuencia procesada
            if sequence_name and sequence and count == 9:
                file_name = f'output_{file_index}.txt'
                with open(file_name, 'w') as out_file:
                    out_file.write(f">|{sequence_name}_ [gene={gene_name}]\n{sequence}\n")

    except FileNotFoundError:
        print(f"File {file_name} not found.")

    return file_names
