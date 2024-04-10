from typing import List, Union
import numpy.typing as npt
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import re

def base_count(fastafile: str) -> List[int]:
    # 課題 1-1   
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq 
    ATGC_count = [seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")]
    return ATGC_count # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    rev_comp_seq = seq.reverse_complement()
    return str(rev_comp_seq)

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    gc_content = []
    for i in range(0, len(seq)-window+1, step):
        gc_content.append((seq[i:i+window].count("G") + seq[i:i+window].count("C")) / window * 100)
    
    return gc_content

def plot_gc_content(gc_content: Union[npt.NDArray[np.float_], List[float]]) -> None:
    # 課題 1-3
    plt.plot(gc_content)
    plt.xlabel("Window")
    plt.ylabel("GC Content (%)")
    plt.savefig("gc_content.png")
    plt.show()

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    rev_com_seq = seq.reverse_complement()

    subseq = []
    if seq.find(motif) != -1:
        subseq.append("F" + str(seq.find(motif)+1))
    if rev_com_seq.find(motif) != -1:
        subseq.append("R" + str(len(seq) - rev_com_seq.find(motif)))

    return subseq

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq = seq_record.seq
    rev_com_seq = seq.reverse_complement()

    trans_list = []
    seq_list = []

    if len(seq) % 3 == 0:
        seq_list.append(seq.translate())
        seq_list.append(seq[1:-2].translate())
        seq_list.append(seq[2:-1].translate())
        seq_list.append(rev_com_seq.translate())
        seq_list.append(rev_com_seq[1:-2].translate())
        seq_list.append(rev_com_seq[2:-1].translate())
    elif len(seq) % 3 == 1:
        seq_list.append(seq[:-1].translate())
        seq_list.append(seq[1:].translate())
        seq_list.append(seq[2:-2].translate())
        seq_list.append(rev_com_seq[:-1].translate())
        seq_list.append(rev_com_seq[1:].translate())
        seq_list.append(rev_com_seq[2:-2].translate())
    else:
        seq_list.append(seq[:-2].translate())
        seq_list.append(seq[1:-1].translate())
        seq_list.append(seq[2:].translate())
        seq_list.append(rev_com_seq[:-2].translate())
        seq_list.append(rev_com_seq[1:-1].translate())
        seq_list.append(rev_com_seq[2:].translate())
    

    for seq in seq_list:
        trans_list.extend(re.findall(r'M[A-LN-Z]*_|M[A-LN-Z]*', str(seq)))
    
    trans_list = list(set(trans_list))
    
    return trans_list

if __name__ == "__main__":
    filepath = "data/NT_113952.1.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath))
    # 課題 1-3 plot
    plot_gc_content(calc_gc_content(filepath))
    # 課題 1-4
    print(search_motif(filepath, "ATG"))
    # 課題 1-5
    print(translate(filepath))
