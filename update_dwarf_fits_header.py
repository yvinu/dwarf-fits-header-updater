#!/usr/bin/env python3
"""
update_dwarf_fits_header.py

Updates the FITS header of DWARF Mini calibration master frames (dark, flat, bias, light)
based on information encoded in the filename, parent folder structure, published DWARF Mini
hardware specifications, and standard FITS header keywords defined in FITS_STANDARD.txt.

DWARF Mini Specifications:
- Telescope Diameter (APTDIA): 30.0 mm
- Focal Length (FOCALLEN): 150.0 mm (f/5.0)
- Sensor (INSTRUME): Sony IMX662
- Native Pixel Size (XPIXSZ, YPIXSZ): 2.90 um
- Plate Scale (SCALE): ~3.987789 arcsec/pixel (at bin 1)
- Cameras: TELE (cam_0) / WIDE (cam_1)
- Filters: Astro (ir_1) / Dual-Band (ir_2)
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

try:
    from astropy.io import fits
except ImportError:
    print("Error: astropy package is required. Install via: pip install astropy", file=sys.stderr)
    sys.exit(1)


# DWARF Mini Hardware Specifications
DWARF_MINI_SPECS = {
    "TELESCOP": "DWARF Mini",
    "APTDIA": 30.0,       # Telescope objective diameter in mm
    "FOCALLEN": 150.0,    # Focal length in mm
    "INSTRUME": "Sony IMX662",
    "NATIVE_PIXEL_SIZE": 2.90,  # um
}

IMAGETYP_MAP = {
    "dark": ("Dark Frame", "Dark"),
    "flat": ("Flat Frame", "Flat"),
    "bias": ("Bias Frame", "Bias"),
    "light": ("Light Frame", "Light"),
}


def parse_filename_and_path(filepath: str | Path) -> Dict[str, Any]:
    """
    Parse metadata encoded in the FITS filename and parent directory path structure.
    Extracts image type, exposure, gain, binning, temperature, stack count, camera module, and filter mode.
    """
    path = Path(filepath)
    filename = path.name
    path_str = str(path)
    parts_lower = [p.lower() for p in path.parts]

    is_cali_frame = "cali_frame" in parts_lower

    # 1. Image Type (dark, flat, bias, light)
    imtype_raw = None
    for cand in ("bias", "flat", "dark", "light"):
        if re.search(r"(?:^|[_\-/\\])" + cand + r"(?:[_\-/\\]|\.fits?)", path_str, re.IGNORECASE):
            imtype_raw = cand
            break

    # If no explicit image type was found:
    # Files inside CALI_FRAME default to 'dark', files outside CALI_FRAME default to 'light'
    if imtype_raw is None:
        imtype_raw = "dark" if is_cali_frame else "light"

    imtype, frame = IMAGETYP_MAP.get(imtype_raw, ("Light Frame", "Light"))

    # 2. Exposure time (EXPTIME)
    exp_m = re.search(r"exp[_\-](?P<exp>[0-9.]+)", filename, re.IGNORECASE) or re.search(r"exp[_\-](?P<exp>[0-9.]+)", path_str, re.IGNORECASE)
    if exp_m:
        exp: Optional[float] = float(exp_m.group("exp"))
    elif imtype_raw == "bias":
        exp = 0.0
    else:
        exp = None

    # 3. Gain (GAIN)
    gain_m = re.search(r"gain[_\-](?P<gain>[0-9.]+)", filename, re.IGNORECASE) or re.search(r"gain[_\-](?P<gain>[0-9.]+)", path_str, re.IGNORECASE)
    if gain_m:
        gain_val = float(gain_m.group("gain"))
        gain: Optional[Union[int, float]] = int(gain_val) if gain_val.is_integer() else gain_val
    else:
        gain = None

    # 4. Binning (XBINNING / YBINNING)
    bin_m = re.search(r"bin[_\-](?P<bin>[0-9]+)", filename, re.IGNORECASE) or re.search(r"bin[_\-](?P<bin>[0-9]+)", path_str, re.IGNORECASE)
    binning = int(bin_m.group("bin")) if bin_m else 1

    # 5. Temperature (CCD-TEMP)
    # Require preceded by non-letter, e.g. 33C, -10C, 33.5C
    # If no temperature is explicitly in filename or folder name, temp is None (do NOT set to 0)
    temp_m = re.search(r"(?<![A-Za-z])(?P<temp>[-+]?[0-9.]+)C(?![A-Za-z])", filename) or re.search(r"(?<![A-Za-z])(?P<temp>[-+]?[0-9.]+)C(?![A-Za-z])", path_str)
    temp: Optional[float] = float(temp_m.group("temp")) if temp_m else None

    # 6. Stack Count (STACKCNT)
    stack_m = re.search(r"stack[_\-](?P<stack>[0-9]+)", filename, re.IGNORECASE) or re.search(r"stack[_\-](?P<stack>[0-9]+)", path_str, re.IGNORECASE)
    stack = int(stack_m.group("stack")) if stack_m else 1

    # 7. Camera (cam_0 ≈ TELE | cam_1 ≈ WIDE)
    cam_name = None
    if re.search(r"cam[_\-]?1|\bwide\b", path_str, re.IGNORECASE):
        cam_name = "WIDE"
    elif re.search(r"cam[_\-]?0|\btele\b", path_str, re.IGNORECASE):
        cam_name = "TELE"
    else:
        cam_name = "TELE"

    # 8. Filter (ir_1 ≈ Astro, ir_2 ≈ Dual-Band)
    # Note: DWARF Mini does not have ir_0
    filter_name = None
    if re.search(r"ir[_\-]?1|\bastro\b", path_str, re.IGNORECASE):
        filter_name = "Astro"
    elif re.search(r"ir[_\-]?2|\bdual\b|\bduo\b|\bband\b|\bnarrow\b", path_str, re.IGNORECASE):
        filter_name = "Dual-Band"

    return {
        "imtype_raw": imtype_raw,
        "IMAGETYP": imtype,
        "FRAME": frame,
        "EXPTIME": exp,
        "GAIN": gain,
        "XBINNING": binning,
        "YBINNING": binning,
        "CCD-TEMP": temp,
        "STACKCNT": stack,
        "CAMNAME": cam_name,
        "FILTER": filter_name,
    }


def parse_filename(filepath: str | Path) -> Dict[str, Any]:
    """Backward compatibility wrapper around parse_filename_and_path."""
    return parse_filename_and_path(filepath)


def calculate_derived_fields(parsed: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate derived fields such as total exposure times, pixel sizes, and scale.
    """
    binning = parsed.get("XBINNING", 1)
    exp = parsed.get("EXPTIME")
    stack = parsed.get("STACKCNT", 1)

    livetime = (exp * stack) if (exp is not None and stack is not None) else None
    # DARKTIME is relevant ONLY for Dark frames
    darktime = livetime if (parsed.get("imtype_raw") == "dark") else None
    pix_sz = DWARF_MINI_SPECS["NATIVE_PIXEL_SIZE"] * binning
    focal_len = DWARF_MINI_SPECS["FOCALLEN"]

    # Plate scale formula: scale (arcsec/pixel) = 206.264806 * pixel_size (um) / focal_length (mm)
    scale = (206.264806 * pix_sz) / focal_len

    return {
        "LIVETIME": livetime,
        "DARKTIME": darktime,
        "XPIXSZ": pix_sz,
        "YPIXSZ": pix_sz,
        "SCALE": scale,
    }


def load_disabled_standard_keys(standard_path: str | Path = "FITS_STANDARD.txt") -> set[str]:
    """
    Parse FITS_STANDARD.txt to find variables prefixed with '/' (which should not be applied).
    Returns a set of uppercase disabled header key names.
    """
    disabled_keys = set()
    p = Path(standard_path)
    if not p.is_file():
        return {"EXPSTART", "EXPEND", "OBJECT", "WB_R", "WB_B", "MIPS-FLO"}

    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if line_str.startswith("/"):
                content = line_str[1:].lstrip()
                if "=" in content:
                    key = content.split("=")[0].strip().upper()
                    if key and re.match(r"^[A-Z0-9_-]{1,8}$", key):
                        disabled_keys.add(key)

    return disabled_keys


def update_fits_header(
    filepath: str | Path,
    output_path: Optional[str | Path] = None,
    standard_path: str | Path = "FITS_STANDARD.txt",
    verbose: bool = True
) -> Path:
    """
    Update the FITS header of the given file to conform to FITS_STANDARD.txt and DWARF Mini specs.
    Variables prefixed with '/' in FITS_STANDARD.txt are excluded and removed.
    """
    filepath = Path(filepath)
    if not filepath.is_file():
        raise FileNotFoundError(f"File not found: {filepath}")

    parsed = parse_filename_and_path(filepath)
    derived = calculate_derived_fields(parsed)
    disabled_keys = load_disabled_standard_keys(standard_path)

    save_path = Path(output_path) if output_path else filepath

    with fits.open(filepath, mode="update" if save_path == filepath else "readonly") as hdul:
        header = hdul[0].header
        data = hdul[0].data

        # Determine dimensions and bitpix
        naxis1 = data.shape[1] if data is not None and data.ndim >= 2 else header.get("NAXIS1", 1920)
        naxis2 = data.shape[0] if data is not None and data.ndim >= 2 else header.get("NAXIS2", 1080)
        
        # Bayer pattern check
        bayerpat = header.get("BAYERPAT", "RGGB    ")

        # BZERO & BSCALE
        bzero = header.get("BZERO", 32768.0 if data is not None and data.dtype.kind == "u" else 0.0)
        bscale = header.get("BSCALE", 1.0)

        now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")

        # Purge DARKTIME from non-dark frames
        if parsed["imtype_raw"] != "dark" and "DARKTIME" in header:
            del header["DARKTIME"]

        # Auto-correct misclassified image type headers on Light frames (e.g. from previous script runs)
        if parsed["imtype_raw"] == "light":
            if header.get("IMAGETYP") in ("Dark Frame", "Dark"):
                header["IMAGETYP"] = (parsed["IMAGETYP"], "Type of image")
            if header.get("FRAME") in ("Dark Frame", "Dark"):
                header["FRAME"] = (parsed["FRAME"], "Frame Type")
            if header.get("OBJECT") in ("Dark Frame", "Dark"):
                header["OBJECT"] = (parsed["FRAME"], "Name of the object of interest")

        # Key-Value pairs with comments matching FITS_STANDARD.txt
        header_updates: list[tuple[str, Any, str]] = [
            ("SIMPLE", True, "file does conform to FITS standard"),
            ("BITPIX", header.get("BITPIX", 16), "number of bits per data pixel"),
            ("NAXIS", 2, "number of data axes"),
            ("NAXIS1", naxis1, "length of data axis 1"),
            ("NAXIS2", naxis2, "length of data axis 2"),
            ("EXTEND", True, "FITS dataset may contain extensions"),
            ("BZERO", float(bzero), "Offset data range to that of unsigned short"),
            ("BSCALE", float(bscale), "Default scaling factor"),
            ("PROGRAM", "DWARF Header Updater v1.2", "Software that created this HDU"),
            ("DATE", now_utc, "UTC date that FITS file was created"),
            ("IMAGETYP", parsed["IMAGETYP"], "Type of image"),
            ("ROWORDER", "TOP-DOWN", "Order of the rows in image array"),
            ("TELESCOP", DWARF_MINI_SPECS["TELESCOP"], "Telescope used to acquire this image"),
            ("FOCALLEN", float(DWARF_MINI_SPECS["FOCALLEN"]), "[mm] Focal length"),
            ("XBINNING", parsed["XBINNING"], "Camera binning mode"),
            ("YBINNING", parsed["YBINNING"], "Camera binning mode"),
            ("XPIXSZ", float(derived["XPIXSZ"]), "[um] Pixel X axis size"),
            ("YPIXSZ", float(derived["YPIXSZ"]), "[um] Pixel Y axis size"),
            ("INSTRUME", DWARF_MINI_SPECS["INSTRUME"], "Instrument name"),
            ("OFFSET", header.get("OFFSET", 0), "Sensor gain offset"),
            ("BAYERPAT", bayerpat, "Bayer color pattern"),
            ("XBAYROFF", 0, "X offset of Bayer array"),
            ("YBAYROFF", 0, "Y offset of Bayer array"),
            ("FRAME", parsed["FRAME"], "Frame Type"),
            ("OBJECT", parsed["FRAME"], "Name of the object of interest"),
            ("APTDIA", float(DWARF_MINI_SPECS["APTDIA"]), "Telescope diameter (mm)"),
            ("SCALE", float(derived["SCALE"]), "arcsecs per pixel"),
        ]

        if parsed["EXPTIME"] is not None:
            header_updates.append(("EXPTIME", float(parsed["EXPTIME"]), "[s] Exposure time duration"))

        if parsed["GAIN"] is not None:
            header_updates.append(("GAIN", parsed["GAIN"], "Sensor gain"))

        # Temperature: Only add if parsed or pre-existing in header
        ccd_temp_val = parsed["CCD-TEMP"] if parsed["CCD-TEMP"] is not None else header.get("CCD-TEMP")
        if ccd_temp_val is not None:
            header_updates.append(("CCD-TEMP", float(ccd_temp_val), "[degC] CCD temperature"))

        # FOCPOS: Only keep if present; do not generate if missing
        if "FOCPOS" in header:
            header_updates.append(("FOCPOS", header["FOCPOS"], "[step] Focuser position"))

        # FOCTEMP: Only add if pre-existing or derived from CCD-TEMP
        if "FOCTEMP" in header:
            header_updates.append(("FOCTEMP", float(header["FOCTEMP"]), "[degC] Focuser temp"))
        elif ccd_temp_val is not None:
            header_updates.append(("FOCTEMP", float(ccd_temp_val), "[degC] Focuser temp (assumed CCD-TEMP)"))

        if parsed["STACKCNT"] is not None:
            header_updates.append(("STACKCNT", parsed["STACKCNT"], "Stack frames"))

        if derived["LIVETIME"] is not None:
            header_updates.append(("LIVETIME", float(derived["LIVETIME"]), "[s] Exposure time after deadtime correction"))

        # DARKTIME: Relevant ONLY for Dark frames
        if derived["DARKTIME"] is not None and parsed["imtype_raw"] == "dark":
            header_updates.append(("DARKTIME", float(derived["DARKTIME"]), "Total Dark Exposure Time (s)"))

        if parsed.get("CAMNAME"):
            header_updates.append(("CAMNAME", parsed["CAMNAME"], "Camera module (TELE = cam_0, WIDE = cam_1)"))

        if parsed.get("FILTER"):
            header_updates.append(("FILTER", parsed["FILTER"], "Filter name"))

        # Filter out disabled keys (variables starting with '/' in FITS_STANDARD.txt)
        header_updates = [item for item in header_updates if item[0] not in disabled_keys]

        # Non-destructive card update:
        # Populate missing standard fields while preserving all pre-existing fields and values intact.
        for key, val, comment in header_updates:
            if key not in header or header[key] is None or str(header[key]).strip() == "":
                header[key] = (val, comment)

        # Standard Comments
        if "FITS (Flexible Image Transport System)" not in str(header.get("COMMENT", "")):
            header.add_comment("FITS (Flexible Image Transport System) format is defined in 'Astronomy")
            header.add_comment("and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")

        header.add_history(f"Updated header to FITS_STANDARD.txt for DWARF Mini on {now_utc}")

        if save_path == filepath:
            hdul.flush()
            if verbose:
                print(f"[SUCCESS] Updated header in-place: {filepath}")
        else:
            hdul.writeto(save_path, overwrite=True)
            if verbose:
                print(f"[SUCCESS] Saved updated FITS file to: {save_path}")

    return save_path


def collect_fits_files(paths: list[str], frame_type: str = "all") -> list[Path]:
    """
    Collect FITS files from a list of paths (files, directories, or glob patterns).
    When traversing directories, filters for the specified frame_type ('dark', 'flat', 'bias', 'light', 'all').
    """
    collected: list[Path] = []
    seen: set[Path] = set()

    def matches_frame_type(p: Path) -> bool:
        if frame_type == "all":
            return True
        meta = parse_filename_and_path(p)
        if meta["imtype_raw"] == frame_type.lower():
            return True
        return False

    for arg in paths:
        p = Path(arg)
        if p.is_dir():
            # Recursively scan directory
            for ext in ("*.fits", "*.fit", "*.FITS", "*.FIT"):
                for f in p.rglob(ext):
                    f_abs = f.resolve()
                    if f_abs not in seen and matches_frame_type(f):
                        seen.add(f_abs)
                        collected.append(f)
        elif p.is_file():
            f_abs = p.resolve()
            if f_abs not in seen:
                seen.add(f_abs)
                collected.append(p)
        else:
            # Handle glob patterns
            parent = p.parent if p.parent.exists() else Path(".")
            matched = list(parent.glob(p.name))
            for f in matched:
                if f.is_file():
                    f_abs = f.resolve()
                    if f_abs not in seen:
                        seen.add(f_abs)
                        collected.append(f)

    return sorted(collected)


def main() -> None:
    default_input = ["CALI_FRAME"] if Path("CALI_FRAME").is_dir() else []

    parser = argparse.ArgumentParser(
        description="Update FITS header of DWARF Mini calibration master frames based on filename, folder structure, and hardware specs."
    )
    parser.add_argument(
        "fits_files",
        nargs="*" if default_input else "+",
        default=default_input,
        help="Path to FITS file(s), directory (e.g., CALI_FRAME), or wildcard pattern (default: CALI_FRAME)"
    )
    parser.add_argument(
        "--type",
        choices=["dark", "flat", "bias", "light", "all"],
        default="all",
        help="Filter frame type when scanning directories (default: all)"
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Modify the FITS file(s) directly in place"
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output filepath (only applicable when processing a single FITS file)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Display header updates without saving changes to disk"
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress detailed console output"
    )

    args = parser.parse_args()
    verbose = not args.quiet

    files_to_process = collect_fits_files(args.fits_files, frame_type=args.type)

    if not files_to_process:
        print(f"[ERROR] No matching FITS files found (frame type: '{args.type}').", file=sys.stderr)
        sys.exit(1)

    if args.output and len(files_to_process) > 1:
        print("[ERROR] --output option can only be used with a single input file.", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Found {len(files_to_process)} '{args.type}' FITS file(s) to process.")

    success_count = 0
    for idx, fits_file in enumerate(files_to_process, 1):
        if verbose:
            print(f"\n[{idx}/{len(files_to_process)}] Processing: {fits_file}")
        try:
            parsed = parse_filename_and_path(fits_file)
            derived = calculate_derived_fields(parsed)
            if verbose:
                print("  Extracted metadata:")
                for k, v in parsed.items():
                    print(f"    - {k}: {v}")
                print("  Calculated derived values:")
                for k, v in derived.items():
                    print(f"    - {k}: {v}")

            if getattr(args, "dry_run", False):
                print("  [DRY RUN] Header would be updated with the above values.")
                success_count += 1
                continue

            in_place = getattr(args, "in_place", False)
            out_path = args.output if args.output else (fits_file if in_place else fits_file.with_name(f"{fits_file.stem}_updated{fits_file.suffix}"))
            update_fits_header(fits_file, output_path=out_path, verbose=verbose)
            success_count += 1

        except Exception as e:
            print(f"[ERROR] Failed to process {fits_file}: {e}", file=sys.stderr)

    if verbose:
        print(f"\nCompleted: Successfully processed {success_count}/{len(files_to_process)} file(s).")


if __name__ == "__main__":
    main()
