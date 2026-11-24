#!/usr/bin/env python3
"""Build a robust stacked star-template spectrum from an object-name catalogue.

Edit the configuration below, then run: python stack_star_template.py
Requires numpy, scipy, astropy, and matplotlib.

FITS extension 8 must contain wavelength, flux, and flux error in rows 0–2.
Input wavelengths must already be in the desired common velocity frame.
The output scatter describes variation among spectra, not propagated error.
"""

import csv
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance
from scipy.interpolate import interp1d

# User configuration
OBJECT_NAME = "M53STAR"  # Exact name in the CSV's Object column.
CATALOG_PATH = Path("object_names.csv")
OUTPUT_DIRECTORY = Path(".")

# Velocity used to position reference lines in the plot (km/s).
# Set the appropriate velocity for the template.
REFERENCE_VELOCITY_KMS = -63.1
REFERENCE_WAVELENGTHS = (5167.3, 5172.7, 5183.6)
PLOT_X_LIMITS = (5100, 5200)  # Use None to display the full wavelength range.
PLOT_Y_LIMITS = (-0.4, 2)  # Use None for automatic flux limits.


def main():
    full_paths = []
    with CATALOG_PATH.open(newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        required = {"Object", "Filename", "Directory"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{CATALOG_PATH} must contain {sorted(required)}.")
        for row in reader:
            if row["Object"] == OBJECT_NAME:
                full_paths.append(os.path.join(row["Directory"], row["Filename"]))

    if not full_paths:
        raise ValueError(f"No files for {OBJECT_NAME!r} found in {CATALOG_PATH}.")

    wavelengths, fluxes = [], []
    for full_path in full_paths:
        with fits.open(full_path) as hdul:
            frame_data = hdul[8].data
            wavelength = np.array(frame_data[0, :], dtype=float, copy=True)
            flux = np.array(frame_data[1, :], dtype=float, copy=True)

        if (wavelength.size < 2 or not np.all(np.isfinite(wavelength))
                or not np.all(np.diff(wavelength) > 0)):
            raise ValueError(
                f"{full_path}: wavelengths must be finite and strictly increasing "
                "with at least two samples."
            )
        wavelengths.append(wavelength)
        fluxes.append(flux)

    min_wave = max(w.min() for w in wavelengths)
    max_wave = min(w.max() for w in wavelengths)
    step = np.mean([np.mean(np.diff(w)) for w in wavelengths])
    if min_wave >= max_wave:
        raise ValueError("The input spectra have no overlapping wavelength range.")
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Wavelength arrays must have a positive wavelength step.")

    common_wavelength = np.arange(min_wave, max_wave, step)
    interp_fluxes = np.array([
        interp1d(wl, fl, kind="linear", bounds_error=False, fill_value=np.nan)(
            common_wavelength
        )
        for wl, fl in zip(wavelengths, fluxes)
    ])

    stacked_flux = biweight_location(interp_fluxes, axis=0, ignore_nan=True)
    stacked_scatter = np.sqrt(
        biweight_midvariance(interp_fluxes, axis=0, ignore_nan=True)
    )
    stacked_data = np.column_stack((
        common_wavelength, stacked_flux, stacked_scatter
    ))

    output_stem = re.sub(r"[^a-zA-Z0-9_-]+", "_", OBJECT_NAME).strip("_").lower()
    if not output_stem:
        raise ValueError("OBJECT_NAME must contain a letter, digit, or hyphen.")
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    output_base = OUTPUT_DIRECTORY / f"{output_stem}_stacked"
    csv_path = output_base.with_suffix(".csv")
    ascii_path = output_base.with_suffix(".ascii")
    png_path = output_base.with_suffix(".png")

    np.savetxt(
        csv_path, stacked_data, fmt="%.16g", delimiter=",",
        header="wavelength,stacked_flux,stacked_scatter", comments=""
    )
    # Columns: wavelength, stacked flux, error among the input spectra.
    np.savetxt(ascii_path, stacked_data, fmt="%.16g", delimiter=" ")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.step(common_wavelength, stacked_flux, color="purple", linewidth=0.8,
            where="mid")
    ax.fill_between(
        common_wavelength, stacked_flux - stacked_scatter,
        stacked_flux + stacked_scatter, alpha=0.4, color="pink", linewidth=0.8
    )
    if REFERENCE_VELOCITY_KMS is not None:
        z = REFERENCE_VELOCITY_KMS / 299792.458
        for line_wave in REFERENCE_WAVELENGTHS:
            ax.axvline(line_wave * (1 + z), linewidth=0.8, color="olive",
                       linestyle="--")
    if PLOT_X_LIMITS is not None:
        ax.set_xlim(*PLOT_X_LIMITS)
    if PLOT_Y_LIMITS is not None:
        ax.set_ylim(*PLOT_Y_LIMITS)
    ax.set_xlabel("Observed Wavelength (Å)")
    ax.set_ylabel("Flux")
    ax.set_title(f"Stacked Star Template: {OBJECT_NAME}")
    fig.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Loaded {len(full_paths)} files for {OBJECT_NAME}:")
    for path in full_paths:
        print(path)
    print("\nSaved files:")
    for path in (png_path, csv_path, ascii_path):
        print(path.resolve())


if __name__ == "__main__":
    main()
