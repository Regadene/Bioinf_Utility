import argparse
from bioinf_utility import filter_fastq
import logging

def parse_args_fastq_fitering():
    parser = argparse.ArgumentParser(description="Filter FASTQ files based on quality, GC content, and length")

    # Input
    parser.add_argument(
        "-i",
        "--inname",
        "--input",
        dest="input",
        type=str,
        required=True,
        help="""Name of the input file in fastq format.""",
    )

    # Output
    parser.add_argument(
        "-o",
        "--outname",
        "--output",
        dest="output",
        default=None,
        type=str,
        required=False,
        help="""Name of the output file.
                If not set, file is saved in the new filtered directory with the same as input file name.""",
    )

    # gc_bounds
    parser.add_argument(
        "--gc",
        "--gc_bounds",
        dest="gc_bounds",
        nargs="+",
        type=int,
        default=[0, 100], 
        required=False,
        help="Sequence GC bounds: provide 1 (max) or 2 integers (min max)"
    )
    
    # length_bounds
    parser.add_argument(
        "-l"
        "--length_bounds",
        dest="length_bounds",
        nargs="+",
        type=int,
        default=[0, 2**32],
        help="Sequence length bounds: provide 1 (max) or 2 integers (min max)")
    
    # quality_threshold
    parser.add_argument(
        "-q"
        "--quality",
        "--quality_threshold",
        dest="quality_threshold",
        type=int,
        default=0,
        help="Minimum quality score: provide integer value")
    
    return parser



def main():
    parser = parse_args_fastq_fitering()
    args = parser.parse_args()
    
    try:
        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=args.gc_bounds,
            length_bounds=args.length_bounds,
            quality_threshold=args.quality_threshold,
        )
    except Exception as e:
        logging.error(str(e))

if __name__ == "__main__":
    main()