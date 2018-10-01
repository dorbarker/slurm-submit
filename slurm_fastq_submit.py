#!/usr/bin/env python
import argparse
import itertools
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('template',
                        type=Path,
                        help='Path to SLURM job template file')

    parser.add_argument('fastqs',
                        type=Path,
                        help='Directory containing FASTQ files')

    parser.add_argument('jobs',
                        type=Path,
                        help='Job script directory')

    parser.add_argument('-o', '--output',
                        required=False,
                        type=Path,
                        help='Output path provided to the template File')

    parser.add_argument('--reference',
                        required=False,
                        type=Path,
                        help='Path to reference genome')

    return parser.parse_args()


def group_fastqs(fastqs: List[Path]) -> Dict[str, Dict[str, str]]:
    """Group FASTQs based on common base filename
    For example, if you have 2 FASTQs:
    - reads_1.fastq
    - reads_2.fastq
    The common name would be `reads` and the files would be
    grouped based on that common name.

    Args:
        fastqs: FASTQ file paths
    Returns:
        list of grouped FASTQs grouped by common base filename
    """

    fastq_regex = re.compile(r'^(.+)\.(fastq|fq)(\.gz|\.bz2)?$')
    illumina_regex  = re.compile(r'_S.+_L\d*1_R\d+_001')

    genome_fastqs = defaultdict(list)

    fwd_rev = {}

    for fastq in fastqs:

        genome_name = re.sub(pattern=illumina_regex,
                              repl='',
                              string=fastq_regex.sub(r'\1', fastq.stem))

        genome_fastqs[genome_name].append(fastq)

    for name, paths in genome_fastqs.items():

        fwd, rev, *_ = sorted(paths)

        fwd_rev[name] = {'fwd': fwd, 'rev': rev}

    print('Grouped {} FASTQs into {} samples'.format(len(fastqs), len(fwd_rev)))

    return fwd_rev


def get_fastq_paths(directory: Path) -> List[Path]:

    extensions = ['*{}{}'.format(a, b)
                  for a, b
                  in itertools.product(['.fq', '.fastq'], ['', '.gz', '.bz2'])]

    globs = [directory.glob(ext) for ext in extensions]

    fastq_paths = [p.absolute() for p in itertools.chain.from_iterable(globs)]

    print('Found {} fastqs in {}'.format(len(fastq_paths), directory))

    return fastq_paths


def main():

    args = arguments()

    job_tmpl = args.template.read_text()

    # get all the *.fastq files in a path
    all_fastq_paths = get_fastq_paths(args.fastqs)

    fastq_paths_names = group_fastqs(all_fastq_paths)

    # writing job files
    for sample_name, paths in fastq_paths_names.items():

        job_filepath = args.jobs / 'slurm-{}.job'.format(sample_name)

        job = job_tmpl.format(genome=sample_name,
                              output=args.output,
                              reference=args.reference,
                              fwd=paths['fwd'],
                              rev=paths['rev'])

        job_filepath.write_text(job)

        print('Wrote {} for sample {} with inputs {} {}'.format(job_filepath,
                                                                sample_name,
                                                                paths['fwd'],
                                                                paths['rev']))


if __name__ == '__main__':
    main()
