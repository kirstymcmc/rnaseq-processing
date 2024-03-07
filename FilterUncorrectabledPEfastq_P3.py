import sys
import gzip
from itertools import zip_longest
import argparse
from os.path import basename

def get_input_streams(r1file, r2file):
    if r1file.endswith('.gz'):
        r1handle = gzip.open(r1file, 'rt', encoding='utf-8')
        r2handle = gzip.open(r2file, 'rt', encoding='utf-8')
    else:
        r1handle = open(r1file, 'r')
        r2handle = open(r2file, 'r')
    
    return r1handle, r2handle

def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector fastq outputs")
    parser.add_argument('-1', '--left_reads', dest='leftreads', type=str, help='R1 fastq file')
    parser.add_argument('-2', '--right_reads', dest='rightreads', type=str, help='R2 fastq file')
    parser.add_argument('-s', '--sample_id', dest='id', type=str, help='sample name to write to log file')
    opts = parser.parse_args()

    r1out = open('unfixrm_%s' % basename(opts.leftreads).replace('.gz', ''), 'w')
    r2out = open('unfixrm_%s' % basename(opts.rightreads).replace('.gz', ''), 'w')

    r1_cor_count = 0
    r2_cor_count = 0
    pair_cor_count = 0
    unfix_r1_count = 0
    unfix_r2_count = 0
    unfix_both_count = 0

    r1_stream, r2_stream = get_input_streams(opts.leftreads, opts.rightreads)

    with r1_stream as f1, r2_stream as f2:
        R1 = grouper(f1, 4)
        R2 = grouper(f2, 4)
        counter = 0
        for entry1, entry2 in zip(R1, R2):
            counter += 1
            if counter % 100000 == 0:
                print(f"{counter} reads processed")

            head1, seq1, placeholder1, qual1 = [i.strip() for i in entry1 if i]
            head2, seq2, placeholder2, qual2 = [j.strip() for j in entry2 if j]

            if 'unfixable' in head1 and 'unfixable' not in head2:
                unfix_r1_count += 1
            elif 'unfixable' in head2 and 'unfixable' not in head1:
                unfix_r2_count += 1
            elif 'unfixable' in [head1, head2]:
                unfix_both_count += 1
            else:
                if 'cor' in head1:
                    r1_cor_count += 1
                if 'cor' in head2:
                    r2_cor_count += 1
                if 'cor' in head1 or 'cor' in head2:
                    pair_cor_count += 1
                
                head1 = head1.split(' ')[0]  # Adjusted to split by space and remove the trailing part
                head2 = head2.split(' ')[0]
                r1out.write('%s\n' % '\n'.join([head1, seq1, placeholder1, qual1]))
                r2out.write('%s\n' % '\n'.join([head2, seq2, placeholder2, qual2]))
    
    total_unfixable = unfix_r1_count + unfix_r2_count + unfix_both_count
    total_retained = counter - total_unfixable

    with open(f'rmunfixable_{opts.id}.log', 'w') as unfix_log:
        unfix_log.write(f'total PE reads:{counter}\nremoved PE reads:{total_unfixable}\nretained PE reads:{total_retained}\nR1 corrected:{r1_cor_count}\nR2 corrected:{r2_cor_count}\npairs corrected:{pair_cor_count}\nR1 unfixable:{unfix_r1_count}\nR2 unfixable:{unfix_r2_count}\nboth reads unfixable:{unfix_both_count}\n')

    r1out.close()
    r2out.close()
    print('done')