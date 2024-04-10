from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO


def enumerate_pairs(fastafile: str) -> List[Tuple[int, int]]:
    # 課題 2-1
    pair_list = []
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq

    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] == "A" and seq[j] == "T":
                pair_list.append((i+1, j+1))
            elif seq[i] == "T" and seq[j] == "A":
                pair_list.append((i+1, j+1))
            elif seq[i] == "C" and seq[j] == "G":
                pair_list.append((i+1, j+1))
            elif seq[i] == "G" and seq[j] == "C":
                pair_list.append((i+1, j+1))
    return pair_list

def enumerate_possible_pairs(fastafile: str, min_distance: int=4) -> List[Tuple[int, int]]:
    # 課題 2-2
    pair_list = []
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq

    min_distance = 4

    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] == "A" and seq[j] == "T" and i + min_distance < j:
                pair_list.append((i+1, j+1))
            elif seq[i] == "T" and seq[j] == "A" and i + min_distance < j:
                pair_list.append((i+1, j+1))
            elif seq[i] == "C" and seq[j] == "G" and i + min_distance < j:
                pair_list.append((i+1, j+1))
            elif seq[i] == "G" and seq[j] == "C" and i + min_distance < j:
                pair_list.append((i+1, j+1))
    return pair_list

def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3
    pair_list = []
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    
    count = 0
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if judge_pairs(seq, i, j, min_length):
                count += 1
                for k in range(1, j-i):
                    if judge_pairs(seq, i+k, j-k, min_length):
                        count += 1
                    else:
                        if count >= min_length:
                            pair_list.append((i+1, j+1, count))
                        count = 0
                        break
            

    return pair_list

def judge_pairs(seq: str, i: int, j: int, min_length: int) -> bool:
    if j - i < min_length:
        return False
    else:
        if seq[i] == "A" and seq[j] == "T":
            return True
        elif seq[i] == "T" and seq[j] == "A":
            return True
        elif seq[i] == "C" and seq[j] == "G":
            return True
        elif seq[i] == "G" and seq[j] == "C":
            return True
        else:
            return False


def create_dotbracket_notation(fastafile: str, min_distance: int=4, min_length: int=2) -> str:
    # 課題 2-4
    pair_list = enumerate_continuous_pairs(fastafile, min_distance, min_length)
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    dotbracket = ["." for _ in range(len(seq))]
    for pair in pair_list:
        for i in range(pair[2]):
            dotbracket[pair[0]+i-1] = "("
            dotbracket[pair[1]-i-1] = ")"

    return "".join(dotbracket)

if __name__ == "__main__":
    #filepath = "data/AUCGCCAU.fasta"
    filepath = "data/NM_014495.4.fasta"
    # 課題 2-1
    print(enumerate_pairs(filepath))
    # 課題 2-2
    print(enumerate_possible_pairs(filepath))
    # 課題 2-3
    print(enumerate_continuous_pairs(filepath, 2))
    # 課題 2-4
    print(create_dotbracket_notation(filepath, 2))


