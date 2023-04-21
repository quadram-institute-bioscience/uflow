def read_fasta(path):
    if path.endswith(".gz"):
        import gzip
    
    seqName = None
    seqComment = None
    with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if seqName is not None:
                    yield seqName, seqComment, sequence
                seqName = line[1:].split()[0]
                seqComment = line[1:].split()[1:] if len(line[1:].split()) > 1 else ""
                sequence = ""
                
            else:
                sequence += line.strip()
    yield seqName, seqComment, sequence

if __name__=="__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        for name, comment, sequence in read_fasta(filename):
            print(f">{name} {comment}\n{sequence}")