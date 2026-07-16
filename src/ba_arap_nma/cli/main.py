from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from ba_arap_nma import __version__
from ba_arap_nma.analysis.modal_space import analyze_modal_space
from ba_arap_nma.core import ca_core as core
from ba_arap_nma.core.structure import Structure, load_structure
from ba_arap_nma.engine.generation import GenerationConfig, generate_branch
from ba_arap_nma.output.fullatom import parse_atom_records
from ba_arap_nma.output.two_state import write_best_comparison
from ba_arap_nma.output.writers import write_branch_outputs, write_json, write_tsv


def parse_int_set(text: str) -> list[int]:
    values: list[int] = []
    for part in text.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            start, end = part.split('-', 1)
            values.extend(range(int(start), int(end) + 1))
        else:
            values.append(int(part))
    return sorted(set(values))


def require_pdb_file(path_value: str, label: str) -> Path:
    path = Path(path_value)
    if path.suffix.lower() != '.pdb':
        raise SystemExit(f'{label} must be a .pdb file: {path}')
    if not path.is_file():
        raise SystemExit(f'{label} does not exist or is not a regular file: {path}')
    return path


def prepare_output_directory(path_value: str, overwrite: bool) -> Path:
    path = Path(path_value)
    if path.exists() and any(path.iterdir()):
        if not overwrite:
            raise SystemExit(
                f'Output directory is not empty: {path}. '
                'Choose a new directory or add --overwrite.'
            )
        resolved = path.resolve()
        if resolved == Path(resolved.anchor):
            raise SystemExit('Refusing to overwrite a filesystem root directory')
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def announce(args: argparse.Namespace, text: str) -> None:
    if not getattr(args, 'quiet', False):
        print(text, flush=True)


def common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--modes', default='1-10', help='Mode ranks, e.g. 1-10 or 1,3,5')
    parser.add_argument(
        '--orientations', nargs='+', type=int, choices=(-1, 1), default=[-1, 1],
        help='One or both orientations: --orientations -1 1'
    )
    parser.add_argument('--max-amplitude', type=float, default=20.0, help='Maximum requested C-alpha RMS amplitude in angstroms')
    parser.add_argument('--amplitude-step', type=float, default=0.25, help='Requested amplitude interval in angstroms')
    parser.add_argument('--cutoff', type=float, default=15.0, help='ANM contact cutoff in angstroms')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--overwrite', action='store_true', help='Replace a non-empty output directory')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress messages')


def validate(
    modes: list[int], orientations: list[int], max_amplitude: float,
    amplitude_step: float, analysis_modes: int | None = None,
    summary_k: int | None = None
) -> None:
    if not modes or min(modes) < 1:
        raise SystemExit('Mode ranks must be positive integers')
    if any(value not in (-1, 1) for value in orientations):
        raise SystemExit('Orientations must be -1 or 1')
    if max_amplitude < 0 or amplitude_step <= 0:
        raise SystemExit('Maximum amplitude must be non-negative and amplitude step must be positive')
    steps = max_amplitude / amplitude_step
    if abs(steps - round(steps)) > 1e-8:
        raise SystemExit('Maximum amplitude must be an integer multiple of amplitude step')
    if analysis_modes is not None and analysis_modes < max(modes):
        raise SystemExit('analysis-modes must be at least the largest generated mode rank')
    if summary_k is not None and analysis_modes is not None and not (1 <= summary_k <= analysis_modes):
        raise SystemExit('summary-k must be between 1 and analysis-modes')


def generate(args: argparse.Namespace) -> None:
    args.pdb = require_pdb_file(args.pdb, 'Input')
    modes = parse_int_set(args.modes)
    orientations = list(args.orientations)
    validate(modes, orientations, args.max_amplitude, args.amplitude_step)
    out = prepare_output_directory(args.output, args.overwrite)

    announce(args, f'[BA-ARAP-NMA] Reading {args.pdb}')
    structure = load_structure(args.pdb, args.chain)
    fullatom_records = parse_atom_records(args.pdb, args.chain)
    announce(args, f'[BA-ARAP-NMA] Calculating constrained modes 1-{max(modes)} for {len(structure.records)} residues')
    basis = core.compute_constraint_compatible_modes(
        structure.coords, structure.segment_ids, 'C2', max(modes), args.cutoff
    )
    if basis.vectors.shape[1] < max(modes):
        raise SystemExit(f'Requested mode rank {max(modes)} exceeds the {basis.vectors.shape[1]} available constrained modes')

    summaries = []
    tasks = [(rank, orientation) for rank in modes for orientation in orientations]
    for index, (rank, orientation) in enumerate(tasks, 1):
        sign = '+' if orientation > 0 else '-'
        announce(args, f'[BA-ARAP-NMA] Branch {index}/{len(tasks)}: mode {rank}, orientation {sign}')
        config = GenerationConfig(rank, orientation, args.cutoff, args.max_amplitude, args.amplitude_step)
        result = generate_branch(structure, config, basis.vectors[:, rank - 1], 'ba_arap_nma')
        write_branch_outputs(
            out / result.branch_id, structure.records, result.frames, result.rows,
            result.summary, result.config, fullatom_records
        )
        summaries.append(result.summary)

    write_tsv(out / 'branch_summary.tsv', summaries)
    write_json(out / 'run_metadata.json', {
        'software': 'BA-ARAP-NMA', 'version': __version__, 'command': 'generate',
        'input_pdb': str(args.pdb), 'chain': args.chain, 'modes': modes,
        'orientations': orientations, 'max_amplitude_A': args.max_amplitude,
        'amplitude_step_A': args.amplitude_step, 'cutoff_A': args.cutoff,
    })
    announce(args, f'[BA-ARAP-NMA] Completed. Results: {out}')


def two_state(args: argparse.Namespace) -> None:
    args.input_pdb = require_pdb_file(args.input_pdb, 'Input state')
    args.target_pdb = require_pdb_file(args.target_pdb, 'Paired second state')
    modes = parse_int_set(args.modes)
    orientations = list(args.orientations)
    validate(modes, orientations, args.max_amplitude, args.amplitude_step, args.analysis_modes, args.summary_k)
    out = prepare_output_directory(args.output, args.overwrite)

    announce(args, '[BA-ARAP-NMA] Mapping and aligning the two PDB states')
    pair = core.prepare_pair(
        args.input_pdb, args.target_pdb, args.input_chain, args.target_chain,
        args.min_identity, args.min_coverage
    )
    structure = Structure(pair.start_records, pair.start_coords, pair.segment_ids, pair.start_parse_audit)
    max_analysis = max(args.analysis_modes, max(modes))
    announce(args, f'[BA-ARAP-NMA] Calculating modal-space analysis through K={max_analysis}')
    rows, modal_summary, linear_basis, constrained_basis = analyze_modal_space(
        structure.coords, pair.target_coords_aligned, structure.segment_ids,
        max_analysis, args.cutoff, args.summary_k
    )
    fullatom_records = parse_atom_records(args.input_pdb, args.input_chain)
    write_tsv(out / 'modal_space_analysis.tsv', rows)
    write_tsv(out / 'modal_space_summary.tsv', [modal_summary])
    write_json(out / 'modal_space_summary.json', modal_summary)

    all_summaries = []
    method_best: dict[str, dict] = {}
    methods = [('ba_arap_nma', constrained_basis), ('linear_ca_anm', linear_basis)]
    total = len(methods) * len(modes) * len(orientations)
    task_index = 0
    for method, basis in methods:
        root = out / method
        root.mkdir(exist_ok=True)
        for rank in modes:
            for orientation in orientations:
                task_index += 1
                sign = '+' if orientation > 0 else '-'
                announce(args, f'[BA-ARAP-NMA] Branch {task_index}/{total}: {method}, mode {rank}, orientation {sign}')
                config = GenerationConfig(rank, orientation, args.cutoff, args.max_amplitude, args.amplitude_step)
                result = generate_branch(
                    structure, config, basis.vectors[:, rank - 1], method,
                    pair.target_coords_aligned
                )
                write_branch_outputs(
                    root / result.branch_id, structure.records, result.frames, result.rows,
                    result.summary, result.config,
                    fullatom_records if method == 'ba_arap_nma' else None
                )
                all_summaries.append(result.summary)
                best_frame = int(result.summary['best_target_frame'])
                item = {
                    'method': method, 'mode_rank': rank, 'orientation': orientation,
                    'frame': best_frame,
                    'requested_amplitude_A': result.summary['best_target_requested_amplitude_A'],
                    'target_rmsd_A': result.summary['best_target_rmsd_A'],
                    'transition_coverage': result.summary['best_transition_coverage'],
                    'input_relative_virtual_angle_mae_deg': result.summary['best_input_relative_virtual_angle_mae_deg'],
                    'coords': result.frames[best_frame],
                }
                if method not in method_best or item['target_rmsd_A'] < method_best[method]['target_rmsd_A']:
                    method_best[method] = item
        write_tsv(root / 'branch_summary.tsv', [row for row in all_summaries if row['method'] == method])

    method_best['_target'] = {'coords': pair.target_coords_aligned}
    write_best_comparison(out, structure.records, pair.target_records, method_best)
    write_tsv(out / 'all_branch_summary.tsv', all_summaries)
    write_json(out / 'run_metadata.json', {
        'software': 'BA-ARAP-NMA', 'version': __version__, 'command': 'two-state',
        'input_pdb': str(args.input_pdb), 'target_pdb': str(args.target_pdb),
        'input_chain': args.input_chain, 'target_chain': args.target_chain,
        'aligned_ca_count': pair.aligned_count, 'sequence_identity': pair.sequence_identity,
        'modes': modes, 'orientations': orientations, 'analysis_modes': max_analysis,
        'summary_k': args.summary_k, 'max_amplitude_A': args.max_amplitude,
        'amplitude_step_A': args.amplitude_step, 'cutoff_A': args.cutoff,
        'target_use': 'retrospective_analysis_only',
    })
    announce(args, f'[BA-ARAP-NMA] Completed. Results: {out}')


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog='ba-arap-nma',
        description='BA-ARAP-NMA 1.0: backbone-angle-compatible as-rigid-as-possible normal-mode analysis'
    )
    parser.add_argument('--version', action='version', version=__version__)
    subparsers = parser.add_subparsers(dest='command', required=True)

    single = subparsers.add_parser('generate', help='Generate amplitude-ordered BA-ARAP-NMA structures from one PDB')
    single.add_argument('--pdb', required=True)
    single.add_argument('--chain', default=None)
    common_options(single)
    single.set_defaults(func=generate)

    paired = subparsers.add_parser('two-state', help='Generate structures and compare two experimentally determined PDB states')
    paired.add_argument('--input-pdb', required=True)
    paired.add_argument('--target-pdb', required=True)
    paired.add_argument('--input-chain', required=True)
    paired.add_argument('--target-chain', required=True)
    paired.add_argument('--analysis-modes', type=int, default=50)
    paired.add_argument('--summary-k', type=int, default=10)
    paired.add_argument('--min-identity', type=float, default=0.85)
    paired.add_argument('--min-coverage', type=float, default=0.80)
    common_options(paired)
    paired.set_defaults(func=two_state)

    version = subparsers.add_parser('version', help='Print software version')
    version.set_defaults(func=lambda _: print(__version__))

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
