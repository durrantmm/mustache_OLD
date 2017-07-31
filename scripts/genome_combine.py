# genome_combine.py <outfile> <file1> <file2> <...>

import sys
import fileinput
from os.path import basename


def rename(f, line, index):
    line = line.partition(">")[2]
    new_line = ">" + f.rpartition(".fasta")[0] + "." + str(index) + " " + line
    return new_line


def concatenate(infasta, dir):
    with open(dir, 'w') as fout:
        for file in infasta:
            index = 1

            with fileinput.input(file) as fin:
                for line in fin:
                    if ">" in line:
                        line = rename(str(basename(file)), line, index)
                        index += 1
                    fout.write(line)
                fout.write("\n")


def main():
    fout = sys.argv[1]
    infasta = sys.argv[2:]
    concatenate(infasta, fout)

if __name__ == "__main__":
    main()