# Vienna RNA I/O 

""" 
Wrapper functions for calling and intepreting ViennaRNA shell commands in Python

Author: David S. White 
Contact: dwhite7@wisc.edu
Updated: 2021-08-03
License: MIT

See https://www.tbi.univie.ac.at for full details of ViennaRNA

Inspired by Arnie (https://github.com/DasLab/arnie)
"""

import os 
import subprocess as sp
import random 
import string 
import shutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

# ------------------------------------------
# Initialize Directories for output
# ------------------------------------------
def init_dir():
    """ 
    Make directories for storing figures output from Vienna RNAfold and RNAcofold. 

    Args: 
        should add to allow path option
    """
    if not os.path.isdir('vienna-figures'):
        os.mkdir('vienna-figures')
        os.mkdir('vienna-figures/RNAfold')
        os.mkdir('vienna-figures/RNAcofold')
    elif not os.path.isdir('vienna-figures/RNAfold'):
        os.mkdir('vienna-figures/RNAfold')
    elif not os.path.isdir('vienna-figures/RNAcofold'):
        os.mkdir('vienna-figures/RNAcofold')
    if not os.path.isdir('temp-files'):
        os.mkdir('temp-files')

#def delete_temp_files():
    # delte all temp files. Might be faster doing in bulk then make and deleting for each run

# ------------------------------------------
# Write Sequence File 
# ------------------------------------------
def write_sequence_file(sequence, constraint=None, file_name=None, file_name_length=10):
    """ 
    Convert sequence(s) to a .txt file to be read

    Args: 
        sequence_1 (str or list of str): nucleic acid sequence
        contraint (str): structure constraints
        file_name (str): name of output seq file. If None, random name generated
        file_name_length (int): length of random letters for file name

    Returns: 
        file_name (str): writes a .seq file with "file_name" in FASTA format for each item in sequence
    """
    if file_name is None:
        file_name = os.getcwd()+'/temp-files/temp_file.txt'

    f = open(file_name,'w') 
    for seq in sequence:
        write_line = '> ' + seq + '\n' + seq
        if constraint is not None: 
            write_line = write_line + '\n' + constraint
        write_line = write_line + '\n' + '\n'
        f.write(write_line)
    f.close()

    return file_name

# ------------------------------------------
# Append Command 
# ------------------------------------------
def append_command(command, motif, constraint, dangles, reweight): 
    """ 
    Append common commands to all the vienna RNA functions 

    Args: 
        command (str): command so far to be modified
        motif (): 
        constraint (bool): add constraint (if True, see write_sequence_file)
        dangles ()

    """
    if motif is not None:
        command.append('--motif="%s"' % motif)

    if constraint is not None:
        command.append('-C')
        # command.append('--enforceConstraint') <- caused [] output

    if not dangles:
        command.append('--dangles=0')
        
    if reweight is not None:
        command.append('--commands=%s' % reweight)

    return command
# ------------------------------------------
# Plot Struture
# ------------------------------------------
def plot_structure(image_name, size_inches=[8,6]):
    """
    Plot .ps secondary structure images from RNAcofold or RNAup using matplotlib.pyplot
    
    Args:
        image_name (str): file name to open
        size_inches (float, 2x2): figure size for ouput 
    """ 
    img = mpimg.imread(image_name)
    plt.imshow(img)
    plt.xticks([])
    plt.yticks([])
    fig = plt.gcf()
    fig.set_size_inches(size_inches[0], size_inches[1])
    plt.show()

def check_image_option(sequence, rna_method, show_image=True, save_image=False):
    """ Plot and/or save image [RNAfold, RNAup]
    Args: 
        sequence (list of str): list of strings of sequences to be plotted
        rna_method (str): sets directory for saving image. can be "RNAfold" or "RNAcofold"
        show_image (bool): plot the image using pyplot
        save_image(bool): save image (as .ps).

    """
    for seq in sequence:
        if len(seq)<42: 
            image_file = seq+'_ss.ps'
        else:
            image_file = seq[0:42]+'_ss.ps'

        if show_image: 
            plot_structure(image_file)

        if save_image:
            target_path = 'vienna-figures/'+rna_method
            if not os.path.isdir(target_path):
                init_dir()
            source = os.getcwd()
            destination = source+'/'+target_path
            shutil.move(os.path.join(source, image_file), os.path.join(destination, image_file))
        else:
            os.remove(image_file)
# ------------------------------------------
# RNAfold
# ------------------------------------------
def RNAfold(sequence, T=37, constraint=None, dangles=True, motif=None, reweight=None, show_image=False, save_image=False): 
    """Predict a minimum free energy structure with Vienna RNAfold

    Args:
        seq (str): nucleic acid sequence
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif
        delete_file (Bool): will delete the sequence file needed for shell command (see write_sequence_file)
        delete_image (Bool): will delete the mfe secondary structure image generated

    Returns:
        dictionary with keys:
            mfe_dot (str): minimum free energy (mfe) secondary structure in dot-bracket notation
            mfe_dG (float): mfe of optimal secondary structure (kcal/mol)
            ensemble_dG (float): the free energy of the thermodynamic ensemble (kcal/mol)
            ensemble_frequency (float): the frequency of the mfe structure in the ensemble

    Notes: 
        > does not work at T == 0

    """
    # write command
    if type(sequence) is str: 
        sequence = [sequence]
    
    sequence_file = write_sequence_file(sequence, constraint)
    command = ['RNAfold', '-T', str(T), '-p0', '−−batch']
    comand = append_command(command, motif, constraint, dangles, reweight)

    if not save_image and not show_image:  
        command.append('--noPS')

    # run shell command
    with open(sequence_file) as f:
        process = sp.Popen(command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = process.communicate()

    if show_image or save_image:
        check_image_option(sequence, 'RNAfold', show_image, save_image)
    
    # interpret output (determined empirically)
    out_text = list(filter(None, stdout.decode('utf-8').split('>')))
    rna_fold_output = list()
    for out_line in out_text:
        out_line = list(filter(None, out_line.split('\n')))
        mfe_dot = out_line[2].split(' ')[0]
        mfe_dG = out_line[2].split(' ')[-1][0:-1]
        if mfe_dG[0]=='(':
            mfe_dG = mfe_dG[1:]
        mfe_dG = float(mfe_dG)
        ensemble_dG = float(out_line[3].split(' ')[-2])
        ensemble_frequency = float(list(filter(None, out_line[4].split(' ')))[-1][0:-1])
        rna_fold_output.append({'mfe_dot': mfe_dot, 'mfe_dG':  mfe_dG, 'ensemble_free_energy': ensemble_dG, 'ensemble_frequency': ensemble_frequency})

    return rna_fold_output

# ------------------------------------------
# RNAup
# ------------------------------------------
def RNAup(sequences, T=37, constraint=None, dangles=True, motif=None, reweight=None, save_file=False): 
    """Predict minimum free energy interaction of two RNAs with Vienna RNAfold.

    Args:
        sequences (list of str): list of at least two nucleic acid sequences (type(sequences[0]) = str).
            > Assumes len(sequences[0]) > len(sequences[1])
            > Both sequences need to be writen 5' -> 3'
        T (float): temperature
        constraint (str): structure constraints
        motif (str): argument to vienna motif
        delete_files (Bool): will delete the sequence file and RNAup output files needed for shell command (see write_sequence_file)

    Returns:
        dictionary with keys:
            mfe_dot (str): minimum free energy (mfe) secondary structure in dot-bracket notation
            seq_1_position (int): positions of sequences[0] in mfe secondary structure (kcal/mol)
            seq_2_position (int): positions of sequences[1] in mfe secondary structure (kcal/mol)
            binding_dG (float): Total free energy of binding (kcal/mol)
            duplex_dG (float): Energy from duplex formation (kcal/mol)
            seq_1_dG (float): Opening energy for the longer sequence (sequences[0]) (kcal/mol)
            seq_2_dG (float): Opening energy for the shorter sequence (sequences[1]) (kcal/mol)

        Notes: 
            > binding_dG = duplex_dG + seq_1_dG + seq_2_dG
            > does not work at T == 0

    """
    # write command 
    sequence_file = write_sequence_file(sequences)
    command = ['RNAup', '-T', str(T), '-b', '-o'] # -o no output file
    command = append_command(command, motif, constraint, dangles, reweight)

    # run shell command
    with open(sequence_file) as f:
        process = sp.Popen(command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = process.communicate()

    # interpret output (determined empirically)
    out_text = list(filter(None, stdout.decode('utf-8').split('>')))[1:]
    rna_up_output = list()
    for out_line in out_text:
        out_line = list(filter(None, out_line.split('\n')[1].split(' ')))
        mfe_dot = out_line[0]
        seq_line_1 = out_line[1].split(',')
        seq_1_idx = list([int(seq_line_1[0])-1, int(seq_line_1[1])])
        seq_line_2 = out_line[3].split(',')
        seq_2_idx = list([int(seq_line_2[0])-1, int(seq_line_2[1])])
        binding_dG = float(out_line[4][1:])
        duplex_dG = float(out_line[6])
        seq_1_dG = float(out_line[8])
        seq_2_dG = float(out_line[10][0:-1])
        rna_up_output.append({'mfe_dot': mfe_dot, 'seq_1_idx': seq_1_idx,  
            'seq_2_idx': seq_2_idx, 'binding_dG': binding_dG, 'duplex_dG':duplex_dG, 'seq_1_dG': seq_1_dG, 'seq_2_dG': seq_2_dG})

    return rna_up_output
    
# ------------------------------------------
# RNAcofold
# ------------------------------------------
def RNAcofold(sequences, T=37, constraint=None, dangles=True, motif=None, reweight=None, show_image=False, save_image=False):
    # write to an output file 
    # generates both a dot plot (if -p) and a a picture of the mfe structures (if -p and/or -p0)

    # reformat -> # notation is 'seq1&seq2' e.g., 'CAGGUAAGUAUA&GUCCAUUCAUA'
    if sequences[0].find('&') > -1:
        combined_sequences = sequences
    else: 
        combined_sequences = list()
        for i in range(1, len(sequences)): # sequences[0] == used for all others
            combined_sequences.append(sequences[i]+'&'+sequences[0])

    sequence_file = write_sequence_file(combined_sequences)
    command = ['RNAcofold', '-T', str(T), '-p0']
    comand = append_command(command, motif, constraint, dangles, reweight)
    if not save_image and not show_image:    #<- option to not even call 
        command.append('--noPS')

    # run shell command
    with open(sequence_file) as f:
        process = sp.Popen(command, stdin=f, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = process.communicate()

    if show_image or save_image:
        check_image_option(combined_sequences, 'RNAcofold', show_image, save_image)

    # Parse output (determined empirically)
    rna_cofold_output = list()
    out_text = list(filter(None, stdout.decode('utf-8').split('>')))
    for out_line in out_text: 
        out_line = out_line.split('\n')
        joint_sequence = out_line[0]
        mfe_dot = out_line[2].split(' ')[0]
        mfe_dG = out_line[2].split(' ')[-1]
        mfe_dG =  mfe_dG.replace('(', '')
        mfe_dG =  mfe_dG.replace(')', '')
        mfe_dG = float(mfe_dG)
        ensemble_dG = float(out_line[3].split(' ')[-2])
        mfe_frequency = float(list(filter(None, out_line[4].split(' ')))[6])
        binding_dG = float(list(filter(None, out_line[4].split(' ')))[-1].split('=')[-1])

        rna_cofold_output.append({'sequences': joint_sequence, 'mfe_dot': mfe_dot, 'mfe_dG': mfe_dG, 'ensemble_dG': ensemble_dG, 
            'mfe_frequency': mfe_frequency, 'binding_dG': binding_dG})

    return rna_cofold_output