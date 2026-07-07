import argparse
import os
import subprocess
import shutil
import time
from gfagraphs import Graph
from pathlib import Path

from pgsimplify.simplify_non_branching_paths import compress_non_branching_paths
from pgsimplify.simplify_snp_mnp import iterative_bubble_compression
from pgsimplify.simplify_small_variants import compress_snarls_pipeline
from pgsimplify.utils import save_offsets
from pgsimplify.offsets_in_gfa import pipeline_offsets

def compress_graph(input_gfa, max_len_to_collapse, tmpdir):
    """
    Compress graph in two steps :
    1. Compress non-branching paths
    2. Compress SNP/MNP chains (only when size is inferior to max_len_to collapse for MNP case)
    
    Parameter
    ----------
    input_gfa : str
        Path to the graph to simplify
    max_len_to_collapse :
        Maximal length for a MNP to simplify
    
    Returns
    -------
    dict[int]
        Offsets dictionnary
    int
        Number of nodes in input graph
    """
    # Load and initialize graph
    graph = Graph(input_gfa)
    graph.compute_neighbors()
    graph.sequence_offsets()

    # Print node number
    nb_nodes_begin = len(graph.segments)
    print(f"Initial node number : {nb_nodes_begin}")

    # Linear compression
    offsets = compress_non_branching_paths(graph)
    nb_nodes_non_branching_path_compression = len(graph.segments)
    nb_removed = nb_nodes_begin - nb_nodes_non_branching_path_compression
    print(f"Linear chains compression : removed {nb_removed} nodes form the graph ({round(nb_removed/nb_nodes_begin*100,ndigits=2)}%)")


    # SNP/MNPs compression
    offsets = iterative_bubble_compression(
        graph,
        offsets,
        max_len_to_collapse
    )
    nb_nodes_snp_mnp_compression = len(graph.segments)
    nb_removed = nb_nodes_begin - nb_nodes_snp_mnp_compression
    print(f"Total SNP/MNPs compression : removed {nb_removed} nodes form the graph ({round(nb_removed/nb_nodes_begin*100,ndigits=2)}%)")

    # Temporary saving the graph
    graph.save_graph(str(tmpdir), minimal=True)

    return offsets, nb_nodes_begin


def compute_snarls(tmpdir: Path):
    """
    Calls vg snarls and outputs snarls in json format
    
    Parameter
    ----------
    tmpdir : Path
        path to compressed graph in which we want to detect snarls
    """
    gfa_file = tmpdir / "compressed_graph.gfa"
    vg_file = tmpdir / "graph.vg"
    snarls_file = tmpdir / "graph.snarls"
    json_file = tmpdir / "graph.json"

    # GFA graph file -> VG graph file
    with open(vg_file, "wb") as out:
        subprocess.run(
            ["vg", "convert", "-f", "-W", str(gfa_file)],
            stdout=out,
            check=True,
        )

    # VG graph file -> binary snarls file
    with open(snarls_file, "wb") as out:
        subprocess.run(
            ["vg", "snarls", str(vg_file)],
            stdout=out,
            check=True,
        )

    # binary snarls file -> json snarls file
    with open(json_file, "w") as out:
        subprocess.run(
            ["vg", "view", "-R", str(snarls_file)],
            stdout=out,
            check=True,
            text=True,
        )

def simplify_graph(input_gfa_file, output_dir, max_len_to_collapse, min_variant_size, save_subgraphs, keep_temp):
    """
    Performs simplification piepline
    
    Parameter
    ----------
    input_gfa : str
        Path to the graph to simplify
    ouput_dir : str
        Directory to store simplified graph and auxilary files
    max_len_to_collapse : int
        Maximal length for a MNP to simplify
    min_variant_size : int
        Variants smaller than min_variant_size will be simplified
    save_subgraphs :
        Tells if subgraphs are saved or not
    """
    # Measure executiont time
    start_time = time.perf_counter()

    # Create output_dir if it doesn't exist yet
    os.makedirs(output_dir, exist_ok=True)

    # Create temporary directory to store temporary files
    tmpdir = Path(output_dir) / "temp"
    os.makedirs(tmpdir, exist_ok=True)

    # Compress graph and store it in temporary directory
    gfa_file = tmpdir / "compressed_graph.gfa"
    offsets, nb_nodes_begin = compress_graph(input_gfa_file, max_len_to_collapse, gfa_file)

    # Save offsets
    offset_file = os.path.join(output_dir, f"offsets.txt")
    save_offsets(offsets, str(offset_file))

    # Compute snarls on compressed graph using vg snarls
    compute_snarls(tmpdir)

    # Simplify small variants
    json_file = tmpdir / "graph.json"
    nb_nodes_end = compress_snarls_pipeline(str(gfa_file), str(json_file), output_dir, min_variant_size, save_subgraphs)

    # Supress temporary directory if the option to keep it is not activated
    if not keep_temp:
        shutil.rmtree(tmpdir)

    # Print simplification summary
    removed_percentage = (nb_nodes_begin - nb_nodes_end) / nb_nodes_begin * 100
    print(f"Removed nodes: {removed_percentage:.2f}% ")

    # Print execution time
    end_time = time.perf_counter()
    print(f"Execution time: {end_time - start_time:.2f} s")


def main():
    parser = argparse.ArgumentParser(
        description="Pangenome graph simplification tool"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Command simplify
    parser_simplify = subparsers.add_parser(
        "simplify",
        help="Simplify small variants in a pangenome graph."
    )

    # Required arguments
    parser_simplify.add_argument(
        "input_gfa",
        type=str,
        help="Input GFA file."
    )

    parser_simplify.add_argument(
        "output_dir",
        type=str,
        help="Output directory where the simplified graph and auxiliary files will be written."
    )

    # Optional arguments
    parser_simplify.add_argument(
        "--max-len-to-collapse",
        type=int,
        default=50,
        help="Maximum MNP length to collapse (default: 50)."
    )

    parser_simplify.add_argument(
        "--min-variant-size",
        type=int,
        default=50,
        help="Minimum variant size to keep during simplification (default: 50)."
    )

    parser_simplify.add_argument(
        "--no-subgraphs-save",
        action="store_false",
        dest="save_subgraphs",
        default=True,
        help="Do not save the subgraphs generated during simplification."
    )

    parser_simplify.add_argument(
        "--keep-temporary-files",
        action="store_true",
        dest="keep_temporary_files",
        default=False,
        help="Keep the temporary fiels generating during piepline."
    )

    # Command offsets
    parser_offset = subparsers.add_parser(
        "offsets",
        help="Compute positions from the output of the simplify command."
    )

    parser_offset.add_argument(
        "input_dir",
        type=str,
        help="Directory produced by the simplify command containing main_graph.gfa, subgraphs/, and offsets.txt."
    )

    parser_offset.add_argument(
        "output",
        type=str,
        help="Output GFA file with computed PO tags."
    )

    args = parser.parse_args()


    if args.command == "simplify":
        simplify_graph(
            args.input_gfa,
            args.output_dir,
            args.max_len_to_collapse,
            args.min_variant_size,
            args.save_subgraphs,
            args.keep_temporary_files
        )

    elif args.command == "offsets":
        pipeline_offsets(
            args.input_dir,
            args.output,
        )


   
if __name__ == "__main__":
    main()