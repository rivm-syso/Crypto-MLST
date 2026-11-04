import argparse
from collections import defaultdict

def base_counts_from_pileup(pile_bases, cons_base):
    counts = defaultdict(int)
    del_count = 0
    handled = 0  # how many "coverage" events we've explicitly handled
    i = 0
    while i < len(pile_bases):
        c = pile_bases[i]
        if c == '^':
            i += 2  # skip '^' and following mapping quality
            continue
        elif c == '$':
            i += 1
            continue
        elif c in ['.', ',']:
            counts[cons_base] += 1
            handled += 1
            i += 1
        elif c.upper() in 'ACGTN':
            counts[c.upper()] += 1
            handled += 1
            i += 1
        elif c == '-':  # Deletion compared to reference
            i += 1
            n = ''
            while i < len(pile_bases) and pile_bases[i].isdigit():
                n += pile_bases[i]
                i += 1
            i += int(n) if n else 0  # skip deleted bases
            del_count += 1
            handled += 1
        elif c == '+':  # Insertion, skip inserted bases
            i += 1
            n = ''
            while i < len(pile_bases) and pile_bases[i].isdigit():
                n += pile_bases[i]
                i += 1
            i += int(n) if n else 0  # skip inserted bases
        else:
            i += 1
    return counts, del_count, handled

def main():
    parser = argparse.ArgumentParser(description="Parse pileup for base support info.")
    parser.add_argument("pileup", help="Input pileup file")
    parser.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    args = parser.parse_args()

    fout = open(args.output, "w") if args.output != "-" else None
    header = "chrom\tpos\tcons_base\tcoverage\tA\tC\tG\tT\tDEL\tOTHER\tfraction"
    if fout:
        print(header, file=fout)
    else:
        print(header)

    with open(args.pileup) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 5:
                continue
            chrom, pos, cons_base_pile, cov, pile_bases = cols[:5]
            coverage = int(cov)
            cons_base_from_fa = cons_base_pile

            if cons_base_from_fa not in "ACGT" or coverage == 0:
                counts = {'A':0, 'C':0, 'G':0, 'T':0}
                del_count = 0
                handled = 0
                frac = 0.0
            else:
                counts, del_count, handled = base_counts_from_pileup(pile_bases, cons_base_from_fa)
                for b in 'ACGT':
                    if b not in counts:
                        counts[b] = 0
                support = counts[cons_base_from_fa]
                frac = support / coverage if coverage > 0 else 0.0
            other = coverage - (counts['A'] + counts['C'] + counts['G'] + counts['T'] + del_count)
            row = f"{chrom}\t{pos}\t{cons_base_from_fa}\t{coverage}\t{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t{del_count}\t{other}\t{frac:.4f}"
            if fout:
                print(row, file=fout)
            else:
                print(row)

    if fout:
        fout.close()

if __name__ == "__main__":
    main()
