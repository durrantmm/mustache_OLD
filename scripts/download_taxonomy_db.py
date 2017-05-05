import sys, os
from os.path import dirname, join, abspath
import wget, tarfile

def main(taxondb):
    os.makedirs(taxondb)
    filename = wget.download("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", taxondb)
    untar(filename, True)


def untar(fname, delete=False):
    directory = os.path.dirname(fname)

    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname)
        tar.extractall(directory)
        tar.close()

        if delete:
            os.remove(fname)
    else:
        print("Not a tar.gz file: %s" % fname)


if __name__ == '__main__':
    taxondb = join( dirname(dirname(abspath(os.sys.argv[0]))), 'taxonomy' )

    main(taxondb)