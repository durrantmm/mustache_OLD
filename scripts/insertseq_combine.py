# is_combine.py <outfile> <file1> <file2> <...>

import sys
import fileinput


# Write array into .txt
def concatenate(infasta, dir):
    with open(dir, 'w') as fout:
        for file in infasta:
            with fileinput.input(file) as fin:
                for line in fin:
                    fout.write(line)
                fout.write("\n")


def main():
    fout = sys.argv[1]
    infasta = sys.argv[2:]
    concatenate(infasta, fout)

if __name__ == "__main__":
    main()
