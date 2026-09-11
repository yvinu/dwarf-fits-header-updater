#!/usr/bin/env python3
"""
update_dwarf_fits_header.py

Updates the FITS header of DWARF Mini calibration master frames (e.g., master darks,
flats, biases) based on information encoded in the filename, published DWARF Mini
hardware specifications, and standard FITS header keywords defined in FITS_STANDARD.txt.

DWARF Mini Specifications:
- Telescope Diameter (APTDIA): 30.0 mm
- Focal Length (FOCALLEN): 150.0 mm (f/5.0)
- Sensor (INSTRUME): Sony IMX662
- Native Pixel Size (XPIXSZ, YPIXSZ): 2.90 um
- Plate Scale (SCALE): ~3.987789 arcsec/pixel (at bin 1)

Filename Convention Example:
dark_exp_15.000000_gain_60_bin_1_33C_stack_20.fits
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

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

# Regex to extract parameters from DWARF Mini calibration master filename
# Matches pattern: {type}_exp_{exposure}_gain_{gain}_bin_{bin}_{temp}C_stack_{stack}.fits
FILENAME_REGEX = re.compile(
    r"(?P<imtype>dark|flat|bias|light)"
    r"[_-]exp[_-](?P<exp>[0-9.]+)"
    r"[_-]gain[_-](?P<gain>[0-9.]+)"
    r"[_-]bin[_-](?P<bin>[0-9]+)"
    r"[_-](?P<temp>[-+]?[0-9.]+)C"
    r"[_-]stack[_-](?P<stack>[0-9]+)",
    re.IGNORECASE
)

IMAGETYP_MAP = {
    "dark": ("Dark Frame", "Dark"),
    "flat": ("Flat Frame", "Flat"),
    "bias": ("Bias Frame", "Bias"),
    "light": ("Light Frame", "Light"),
}


def parse_filename(filepath: str | Path) -> Dict[str, Any]:
    """
    Parse metadata encoded in the FITS filename.
    Returns a dictionary of extracted parameters.
    """
    filename = Path(filepath).name
    match = FILENAME_REGEX.search(filename)
    if not match:
        raise ValueError(
            f"Filename '{filename}' does not match expected DWARF calibration pattern: "
            "'{type}_exp_{exp}_gain_{gain}_bin_{bin}_{temp}C_stack_{stack}.fits'"
        )

    imtype_raw = match.group("imtype").lower()
    imtype, frame = IMAGETYP_MAP.get(imtype_raw, ("Dark Frame", "Dark"))

    exp = float(match.group("exp"))
    gain_val = float(match.group("gain"))
    gain = int(gain_val) if gain_val.is_integer() else gain_val
    binning = int(match.group("bin"))
    temp = float(match.group("temp"))
    stack = int(match.group("stack"))

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
    }


def calculate_derived_fields(parsed: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate derived fields such as total exposure times, pixel sizes, and scale.
    """
    binning = parsed["XBINNING"]
    exp = parsed["EXPTIME"]
    stack = parsed["STACKCNT"]

    livetime = exp * stack
    darktime = livetime
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

    parsed = parse_filename(filepath)
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

        # Key-Value pairs with comments matching FITS_STANDARD.txt
        header_updates = [
            ("SIMPLE", True, "file does conform to FITS standard"),
            ("BITPIX", header.get("BITPIX", 16), "number of bits per data pixel"),
            ("NAXIS", 2, "number of data axes"),
            ("NAXIS1", naxis1, "length of data axis 1"),
            ("NAXIS2", naxis2, "length of data axis 2"),
            ("EXTEND", True, "FITS dataset may contain extensions"),
            ("BZERO", float(bzero), "Offset data range to that of unsigned short"),
            ("BSCALE", float(bscale), "Default scaling factor"),
            ("PROGRAM", "DWARF Header Updater v1.0", "Software that created this HDU"),
            ("DATE", now_utc, "UTC date that FITS file was created"),
            ("IMAGETYP", parsed["IMAGETYP"], "Type of image"),
            ("ROWORDER", "TOP-DOWN", "Order of the rows in image array"),
            ("EXPTIME", float(parsed["EXPTIME"]), "[s] Exposure time duration"),
            ("TELESCOP", DWARF_MINI_SPECS["TELESCOP"], "Telescope used to acquire this image"),
            ("FOCALLEN", float(DWARF_MINI_SPECS["FOCALLEN"]), "[mm] Focal length"),
            ("XBINNING", parsed["XBINNING"], "Camera binning mode"),
            ("YBINNING", parsed["YBINNING"], "Camera binning mode"),
            ("XPIXSZ", float(derived["XPIXSZ"]), "[um] Pixel X axis size"),
            ("YPIXSZ", float(derived["YPIXSZ"]), "[um] Pixel Y axis size"),
            ("INSTRUME", DWARF_MINI_SPECS["INSTRUME"], "Instrument name"),
            ("CCD-TEMP", float(parsed["CCD-TEMP"]), "[degC] CCD temperature"),
            ("GAIN", parsed["GAIN"], "Sensor gain"),
            ("OFFSET", header.get("OFFSET", 0), "Sensor gain offset"),
            ("BAYERPAT", bayerpat, "Bayer color pattern"),
            ("XBAYROFF", 0, "X offset of Bayer array"),
            ("YBAYROFF", 0, "Y offset of Bayer array"),
        ]

        # FOCPOS: Only keep if present; do not generate if missing
        if "FOCPOS" in header:
            header_updates.append(("FOCPOS", header["FOCPOS"], "[step] Focuser position"))

        # FOCTEMP: Assumed to be same temperature as sensor (CCD-TEMP) if missing
        foctemp_val = float(header["FOCTEMP"]) if "FOCTEMP" in header else float(parsed["CCD-TEMP"])
        header_updates.append((
            "FOCTEMP",
            foctemp_val,
            "[degC] Focuser temp (assumed CCD-TEMP)"
        ))

        header_updates.extend([
            ("STACKCNT", parsed["STACKCNT"], "Stack frames"),
            ("LIVETIME", float(derived["LIVETIME"]), "[s] Exposure time after deadtime correction"),
            ("OBJECT", parsed["FRAME"], "Name of the object of interest"),
            ("DARKTIME", float(derived["DARKTIME"]), "Total Dark Exposure Time (s)"),
            ("FRAME", parsed["FRAME"], "Frame Type"),
            ("APTDIA", float(DWARF_MINI_SPECS["APTDIA"]), "Telescope diameter (mm)"),
            ("SCALE", float(derived["SCALE"]), "arcsecs per pixel"),
        ])

        # Filter out disabled keys (variables starting with '/' in FITS_STANDARD.txt).
        # Disabled keys are not automatically generated/applied by the script,
        # but if pre-existing in the original file header, they are preserved intact.
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


def collect_fits_files(paths: list[str], frame_type: str = "dark") -> list[Path]:
    """
    Collect FITS files from a list of paths (files, directories, or glob patterns).
    When traversing directories, filters for the specified frame_type (e.g. 'dark').
    """
    collected: list[Path] = []
    seen: set[Path] = set()

    def matches_frame_type(p: Path) -> bool:
        if frame_type == "all":
            return True
        filename_lower = p.name.lower()
        parent_parts = [part.lower() for part in p.parts]
        # Check filename prefix or parent directory names
        if filename_lower.startswith(f"{frame_type}_") or frame_type in parent_parts:
            return True
        # Try regex parse
        match = FILENAME_REGEX.search(p.name)
        if match and match.group("imtype").lower() == frame_type:
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
        description="Update FITS header of DWARF Mini calibration master frames based on filename and hardware specs."
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
        default="dark",
        help="Filter frame type when scanning directories (default: dark)"
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
            parsed = parse_filename(fits_file)
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
