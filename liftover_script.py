import gzip

def read_lines(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            yield line.rstrip('\n')

def read_liftover_lines(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            yield line.rstrip('\n')

def write_gz_line(file, line):
    file.write((line + '\n').encode('utf-8'))

def main():
    # File paths
    input_file = "Raw_Data/EA4_additive_excl_23andMe.txt"
    liftover_bed = "Raw_Data/hglft_genome_315c7d_19acb0.bed" #obtained from the liftover page : https://genome.ucsc.edu/cgi-bin/hgLiftOver?token=0.z-mXvZ0uBQt24k5jrd7hkv5PMFlkevt0XwnWnaXPkZZp5t9mDQmHoINu15ch33NB1UFj4kAIdHHni9PooxjVNUWlizInW24tAEf7borcB-iOY2qJU2ZlqCc3AW-Og4BXcZNsz7MquWUKVNJCwk1u5lQZ2dLt0xCZ04vLBq9yvlltXayTBKPo-r2KybsP5iZqTSHonz1o4rIs0MgYWiTipwD8ahxvF8WlRTU5QHJ5AfdxkUn4TmxQFgX0213krBA4TRa6-nGC2Y_nKiF6g_6LbTCH-sOwsdYS6Sqv40itBszpDofDJeDmGgucXA-DsuZ7yMC8dT-sb5X3BVGE0qA5vkfUmNI61T-y5X6toGSvTz3IR6Ir8gggfXUdSzPcIar-fuifexPovTLLUp84g0HEzYKOqN4RX4Ti-EtS06IHlhl8tpyRgFBs6iUpQ399gwmf6KonWoILxXKfftA8b6z7a6fGYDbuPoycEL8HNeFrhSg_SXNvfcSfLKGH_NqmaYcmmjq3-0zinlZVRvL4EYIWiADtQncFXfIHhm1dJzKLT_4cjExuWIqNR1q6b43PBEj0FlV7g5ofpZKgzEIdc6XScavWcqRHxynjkHVndcLI1zeatoOkHOMC-cPJCE07z5m7xBxwdJ4C1rFM2W1zZJ3lW1Em-sMJ4ughFCqRVmHJEs9GKy_Qwd-fBy1JaFHPGQ9B1ZACAOepBlV9NqVJTlD6FtSck1PnBigbEjj8P2pZZQGDrOgcFuf5JidWmgPLUC4Ou0tP0QoR3FK2cQcpfZAE-fIN5hs9XMigbvKza1laYGzKHTI-gYj6ZT6tBuo20bbcnGUB_ytUPc1aIM2zZcpK0ZOASyhZ7_QcMkbN9hvDlto1iWgpxQU5K041unOpWUIg.hx7nML5Kfm2EOZkmhL4RLw.eae13f0a0de77765d9cc1905d36994a23f19993fda2f2ae634f2efddb7210ab9
    output_file = "Raw_Data/EA4_additive_excl_23andMe_liftover.tsv.gz"

    rf_lines = read_lines(input_file)
    lf_lines = read_liftover_lines(liftover_bed)

    # Read header
    row_in_rf = next(rf_lines)
    rf_iter = iter(rf_lines)
    lf_iter = iter(lf_lines)

    # Open output file
    with gzip.open(output_file, 'wb') as wf:
        # Write header
        write_gz_line(wf, row_in_rf)

        c = -1
        for liftover_line in lf_iter:
            split = liftover_line.split('\t')
            r = int(split[3])
            if r % 1000 == 0:
                print(r, file=sys.stderr)

            while c != r:
                try:
                    row_in_rf = next(rf_iter)
                except StopIteration:
                    print("Reached end of input before expected.")
                    return
                c += 1

            # Modify the position
            splito = row_in_rf.split('\t')
            splito[2] = split[2]  # New position from liftover BED

            # Reconstruct line
            new_line = '\t'.join(splito)
            write_gz_line(wf, new_line)

if __name__ == "__main__":
    import sys
    main()
