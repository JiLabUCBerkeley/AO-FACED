"""
Register a 3-D TIFF stack with CaImAn motion correction.

The input stack is treated as a movie with shape (frames, y, x). For an AO
volume, this means each z slice is registered in XY. By default, a max
projection is prepended as frame 0 so CaImAn has a stable reference, then that
reference frame is removed from the saved registered stack.

Example:
    python caiman_register_tiff_stack.py
    python caiman_register_tiff_stack.py Img_sampleAO_avg_noInterCorr.tif Img_sysAO_avg_noInterCorr.tif
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np
import tifffile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Register a 3-D TIFF stack using CaImAn MotionCorrect."
    )
    parser.add_argument(
        "input_tifs",
        nargs="*",
        default=[
            "Img_sampleAO_avg_noInterCorr.tif",
            "Img_sysAO_avg_noInterCorr.tif",
        ],
        help=(
            "Input 3-D TIFF stack(s). Default: Img_sampleAO_avg_noInterCorr.tif "
            "and Img_sysAO_avg_noInterCorr.tif"
        ),
    )
    parser.add_argument(
        "--output-tif",
        default=None,
        help="Registered TIFF output. Only valid with one input. Default: <input stem>_caiman_reg.tif",
    )
    parser.add_argument(
        "--shifts-csv",
        default=None,
        help="CSV output for CaImAn shifts. Only valid with one input. Default: <input stem>_caiman_shifts.csv",
    )
    parser.add_argument(
        "--no-reference-frame",
        action="store_true",
        help="Do not prepend a max-projection reference frame before correction.",
    )
    parser.add_argument(
        "--pw-rigid",
        action="store_true",
        help="Use piecewise-rigid registration instead of global rigid shifts.",
    )
    parser.add_argument(
        "--max-shift-yx",
        type=int,
        nargs=2,
        default=(15, 15),
        metavar=("Y", "X"),
        help="Maximum allowed rigid shift in pixels as Y X. Default: 15 15",
    )
    parser.add_argument(
        "--strides-yx",
        type=int,
        nargs=2,
        default=(96, 96),
        metavar=("Y", "X"),
        help="Patch stride in pixels for piecewise-rigid mode. Default: 96 96",
    )
    parser.add_argument(
        "--overlaps-yx",
        type=int,
        nargs=2,
        default=(48, 48),
        metavar=("Y", "X"),
        help="Patch overlap in pixels for piecewise-rigid mode. Default: 48 48",
    )
    parser.add_argument(
        "--max-deviation-rigid",
        type=int,
        default=5,
        help="Max patch deviation from rigid shift for piecewise mode. Default: 5",
    )
    parser.add_argument(
        "--n-processes",
        type=int,
        default=1,
        help="Number of local worker processes for CaImAn. Default: 1",
    )
    return parser.parse_args()


def make_reference_augmented_tif(input_tif: Path, stack: np.ndarray) -> Path:
    reference = np.max(stack, axis=0, keepdims=True)
    augmented = np.concatenate([reference, stack], axis=0)
    augmented_tif = input_tif.with_name(f"{input_tif.stem}_with_maxref_tmp_{os.getpid()}.tif")
    tifffile.imwrite(augmented_tif, augmented)
    return augmented_tif


def load_caiman_memmap(memmap_file: str | Path) -> np.ndarray:
    import caiman as cm

    yr, dims, t_frames = cm.load_memmap(str(memmap_file))
    return np.reshape(yr.T, (t_frames, *dims), order="F")


def write_rigid_shift_csv(csv_path: Path, shifts: np.ndarray, dropped_reference: bool) -> None:
    shifts = np.asarray(shifts)
    if dropped_reference:
        shifts = shifts[1:]

    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frame_index", "shift_y_px", "shift_x_px"])
        for frame_index, shift in enumerate(shifts):
            writer.writerow([frame_index, float(shift[0]), float(shift[1])])


def write_piecewise_shift_csv(
    csv_path: Path,
    y_shifts: np.ndarray,
    x_shifts: np.ndarray,
    dropped_reference: bool,
) -> None:
    y_shifts = np.asarray(y_shifts)
    x_shifts = np.asarray(x_shifts)
    if dropped_reference:
        y_shifts = y_shifts[1:]
        x_shifts = x_shifts[1:]

    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "frame_index",
                "mean_shift_y_px",
                "mean_shift_x_px",
                "patch_shifts_y_px_json",
                "patch_shifts_x_px_json",
            ]
        )
        for frame_index, (shift_y, shift_x) in enumerate(zip(y_shifts, x_shifts)):
            writer.writerow(
                [
                    frame_index,
                    float(np.mean(shift_y)),
                    float(np.mean(shift_x)),
                    json.dumps(np.asarray(shift_y, dtype=float).tolist()),
                    json.dumps(np.asarray(shift_x, dtype=float).tolist()),
                ]
            )


def register_stack(
    args: argparse.Namespace,
    input_tif: Path,
    output_tif: Path,
    shifts_csv: Path,
    dview: object | None,
) -> None:
    stack = tifffile.imread(input_tif)
    if stack.ndim != 3:
        raise ValueError(f"Expected a 3-D stack with shape (frames, y, x), got {stack.shape}")

    dropped_reference = not args.no_reference_frame
    caiman_input = input_tif
    if dropped_reference:
        caiman_input = make_reference_augmented_tif(input_tif, stack)

    from caiman.motion_correction import MotionCorrect
    from caiman.source_extraction.cnmf import params as params

    opts_dict = {
        "fnames": [str(caiman_input)],
        "is3D": False,
        "pw_rigid": args.pw_rigid,
        "max_shifts": tuple(args.max_shift_yx),
        "strides": tuple(args.strides_yx),
        "overlaps": tuple(args.overlaps_yx),
        "max_deviation_rigid": args.max_deviation_rigid,
        "border_nan": "copy",
    }
    opts = params.CNMFParams(params_dict=opts_dict)

    try:
        print(f"Registering: {input_tif}")
        mc = MotionCorrect([str(caiman_input)], dview=dview, **opts.get_group("motion"))
        mc.motion_correct(save_movie=True)

        memmap_file = mc.fname_tot_els[0] if args.pw_rigid else mc.fname_tot_rig[0]
        registered = load_caiman_memmap(memmap_file)

        if dropped_reference:
            registered = registered[1:]

        tifffile.imwrite(output_tif, registered.astype(stack.dtype, copy=False))
        if args.pw_rigid:
            write_piecewise_shift_csv(
                shifts_csv,
                np.asarray(mc.y_shifts_els),
                np.asarray(mc.x_shifts_els),
                dropped_reference,
            )
        else:
            write_rigid_shift_csv(shifts_csv, np.asarray(mc.shifts_rig), dropped_reference)

        print(f"Input stack:      {input_tif}")
        print(f"Registered TIFF: {output_tif}")
        print(f"Shift CSV:       {shifts_csv}")
        print(f"Registered shape:{registered.shape}")
    finally:
        if dropped_reference and caiman_input.exists():
            try:
                caiman_input.unlink()
            except PermissionError:
                print(f"Temporary file still in use, leaving it in place: {caiman_input}")


def main() -> None:
    args = parse_args()
    if len(args.input_tifs) > 1 and (args.output_tif or args.shifts_csv):
        raise ValueError("--output-tif and --shifts-csv can only be used with one input TIFF")

    import caiman as cm

    dview = None
    if args.n_processes > 1:
        cluster = cm.cluster.setup_cluster(
            backend="local",
            n_processes=args.n_processes,
            single_thread=False,
        )
        if len(cluster) == 3:
            _, dview, _ = cluster
        else:
            dview = cluster[0]

    try:
        for input_name in args.input_tifs:
            input_tif = Path(input_name).resolve()
            output_tif = Path(
                args.output_tif or input_tif.with_name(f"{input_tif.stem}_caiman_reg.tif")
            )
            shifts_csv = Path(
                args.shifts_csv or input_tif.with_name(f"{input_tif.stem}_caiman_shifts.csv")
            )
            register_stack(args, input_tif, output_tif, shifts_csv, dview)
    finally:
        if dview is not None:
            cm.stop_server(dview=dview)


if __name__ == "__main__":
    main()
