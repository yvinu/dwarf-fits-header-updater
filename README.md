# DWARF Mini FITS Header Updater

A Python utility for standardizing and populating FITS headers of master calibration frames (darks, flats, biases) created by the **DWARF Mini smart telescope** to match astronomical standards and software expectations (such as Siril, PixInsight, and Astro Pixel Processor).

## Features

- **Automated Filename Metadata Parsing**: Extracts `IMAGETYP`, `FRAME`, exposure duration (`EXPTIME`), gain (`GAIN`), binning (`XBINNING`/`YBINNING`), temperature (`CCD-TEMP`), and stack count (`STACKCNT`) from standard DWARF filenames (e.g. `dark_exp_15.000000_gain_60_bin_1_33C_stack_20.fits`).
- **DWARF Mini Spec Alignment**: Automatically populates physical hardware specifications:
  - Objective Aperture (`APTDIA`): **30.0 mm**
  - Focal Length (`FOCALLEN`): **150.0 mm** (f/5.0)
  - Camera Sensor (`INSTRUME`): **Sony IMX662**
  - Native Pixel Size (`XPIXSZ`, `YPIXSZ`): **2.90 µm** ($\text{bin} \times 2.90\,\mu\text{m}$)
  - Plate Scale (`SCALE`): **3.987789 arcsec/pixel** (at bin 1)
- **Derived Computations**: Automatically calculates total integration times (`LIVETIME` & `DARKTIME` = $\text{EXPTIME} \times \text{STACKCNT}$) and optical plate scale.
- **Dynamic Standard Alignment (`FITS_STANDARD.txt`)**: Respects configuration standards defined in `FITS_STANDARD.txt`. Variables prefixed with `/` (e.g. `/EXPSTART`, `/EXPEND`, `/OBJECT`, `/WB_R`, `/WB_B`) are not auto-generated if missing.
- **Non-Destructive Update Policy**: Preserves all pre-existing populated header cards, custom metadata fields, and historical comments without deleting or overwriting them.
- **Recursive Directory Traversal**: Scans calibration folder structures (such as `CALI_FRAME/`) and filters by frame type (`dark`, `flat`, `bias`, `light`, `all`).

---

## Installation

### Prerequisites

- Python 3.8+
- `astropy` library

Install dependencies via pip:

```bash
pip install astropy
```

---

## Usage

### 1. Dry Run (Preview Changes)

Preview extracted parameters and header updates without altering disk files:

```bash
python3 update_dwarf_fits_header.py dark_exp_15.000000_gain_60_bin_1_33C_stack_20.fits --dry-run
```

### 2. Single File Update

Update a FITS file header in-place:

```bash
python3 update_dwarf_fits_header.py dark_exp_15.000000_gain_60_bin_1_33C_stack_20.fits --in-place
```

Or save as a new updated file (`_updated.fits`):

```bash
python3 update_dwarf_fits_header.py dark_exp_15.000000_gain_60_bin_1_33C_stack_20.fits
```

### 3. Calibration Directory Batch Update

Traverse a calibration library structure (e.g., `CALI_FRAME/`) and update all dark master headers in place:

```bash
python3 update_dwarf_fits_header.py CALI_FRAME --in-place
```

### 4. Filter by Frame Type

Filter frame types when scanning directories (`dark`, `flat`, `bias`, `light`, or `all`):

```bash
python3 update_dwarf_fits_header.py CALI_FRAME --type dark --in-place
```

---

## Command-Line Options

```text
usage: update_dwarf_fits_header.py [-h] [--type {dark,flat,bias,light,all}] [--in-place] [--output OUTPUT] [--dry-run] [--quiet] [fits_files ...]

positional arguments:
  fits_files            Path to FITS file(s), directory (e.g., CALI_FRAME), or wildcard pattern (default: CALI_FRAME)

options:
  -h, --help            show this help message and exit
  --type {dark,flat,bias,light,all}
                        Filter frame type when scanning directories (default: dark)
  --in-place            Modify the FITS file(s) directly in place
  --output OUTPUT, -o OUTPUT
                        Output filepath (only applicable when processing a single FITS file)
  --dry-run             Display header updates without saving changes to disk
  --quiet, -q           Suppress detailed console output
```

---

## FITS Header Standard (`FITS_STANDARD.txt`)

The script reads `FITS_STANDARD.txt` to enforce consistent header card names, ordering, and comments.

Variables prefixed with `/` in `FITS_STANDARD.txt` indicate fields that should **not** be automatically added if missing (such as target-specific keywords for dark calibration masters). If a file already contains pre-existing populated values for these keys, they are safely preserved.

---

## License

MIT License.
