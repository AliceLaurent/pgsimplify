import argparse
import os
import subprocess
import shutil
import time
import resource
from gfagraphs import Graph
from pathlib import Path

from pgsimplify.simplify_non_branching_paths import compress_non_branching_paths
from pgsimplify.simplify_snp_mnp import compress_bubbles_chains
from pgsimplify.simplify_small_variants import compress_snarls_pipeline
from pgsimplify.offsets_in_gfa import pipeline_offsets
from pgsimplify.utils import write_report

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
    compress_non_branching_paths(graph)
    nb_nodes_non_branching_path_compression = len(graph.segments)
    nb_removed = nb_nodes_begin - nb_nodes_non_branching_path_compression
    print(f"Non-branching paths compression : removed {nb_removed} nodes form the graph ({round(nb_removed/nb_nodes_begin*100,ndigits=2)}%)")
    nb_nodes_middle = len(graph.segments)

    # SNP/MNPs compression
    compress_bubbles_chains(
                graph=graph,
                max_len=max_len_to_collapse
    )
    nb_nodes_snp_mnp_compression = len(graph.segments)
    nb_removed = nb_nodes_middle - nb_nodes_snp_mnp_compression
    print(f"SNP/MNPs compression : removed {nb_removed} nodes form the graph ({round(nb_removed/nb_nodes_middle*100,ndigits=2)}%)")
    print(f"Number of nodes after non-branching paths and SNP/MNPs compression : {nb_nodes_snp_mnp_compression}")

    # Temporary saving the graph
    graph.save_graph(str(tmpdir), minimal=True)

    return {
    "initial_nodes": nb_nodes_begin,
    "after_non_branching": nb_nodes_non_branching_path_compression,
    "after_snp": nb_nodes_snp_mnp_compression,
    "removed_non_branching": nb_nodes_begin - nb_nodes_non_branching_path_compression,
    "removed_snp": nb_nodes_middle - nb_nodes_snp_mnp_compression,
}


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

    input_gfa_file = Path(input_gfa_file)
    output_dir = Path(output_dir)

    # Vérification du fichier d'entrée
    if not input_gfa_file.is_file():
        raise FileNotFoundError(
            f"Input GFA file not found: {input_gfa_file}"
        )

    if input_gfa_file.suffix.lower() != ".gfa":
        raise ValueError(
            f"Input must be a .gfa file, got: {input_gfa_file}"
        )

    # Vérification du dossier de sortie
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(
            f"Output path exists but is not a directory: {output_dir}"
        )
    
    # Create output_dir if it doesn't exist yet
    os.makedirs(output_dir, exist_ok=True)

    # Create temporary directory to store temporary files
    tmpdir = Path(output_dir) / "temp"
    os.makedirs(tmpdir, exist_ok=True)

    # Compress graph and store it in temporary directory
    gfa_file = tmpdir / "compressed_graph.gfa"
    compression_stats = compress_graph(str(input_gfa_file), max_len_to_collapse, gfa_file)

    # Compute snarls on compressed graph using vg snarls
    compute_snarls(tmpdir)

    # Simplify small variants
    json_file = tmpdir / "graph.json"
    snarl_stats = compress_snarls_pipeline(str(gfa_file), str(json_file), str(output_dir), min_variant_size, save_subgraphs)

    # Supress temporary directory if the option to keep it is not activated
    if not keep_temp:
        shutil.rmtree(tmpdir)

    # Print simplification summary
    removed_percentage = (compression_stats['initial_nodes'] - snarl_stats["after_small_variants"]) / compression_stats['initial_nodes'] * 100
    print(f"Removed nodes: {removed_percentage:.2f}% ")

    # Print execution time
    elapsed = time.perf_counter() - start_time

    # Memory peak
    peak_ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    # ru_maxrss in Ko for Linux
    peak_ram_gb = peak_ram / 1024 / 102

    write_report(
        output_dir=output_dir,
        input_gfa=input_gfa_file,
        max_len=max_len_to_collapse,
        min_variant=min_variant_size,
        save_subgraphs=save_subgraphs,
        compression=compression_stats,
        snarls=snarl_stats,
        elapsed=elapsed,
        peak_ram_gb=peak_ram_gb,
    )

    print(f"Execution time: {elapsed:.2f} s")
    print(f"Peak RAM: {peak_ram_gb:.2f} GB")

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
        help="Directory produced by the simplify command containing main_graph.gfa and subgraphs/ directory containing the subgraphs"
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