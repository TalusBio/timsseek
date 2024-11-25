import numpy as np
from typing import List, Tuple
import sys
from elfragmentadonnx.model import OnnxPeptideTransformer
from typing import List, Tuple, Generator
from tqdm.auto import tqdm
from rich.pretty import pprint
import rustyms
from pyteomics import fasta, parser
import json

# This makes a lot slimmer the logging ingofrmation

model = OnnxPeptideTransformer.default_model()


def convolve(a: List[float], b: List[float]) -> List[float]:
    """
    Performs a custom convolution operation on two arrays of length 4.

    Args:
        a: First array of 4 floating point numbers
        b: Second array of 4 floating point numbers

    Returns:
        List of 4 floating point numbers representing the convolution result
    """
    return [
        a[0] * b[0],
        a[0] * b[1] + a[1] * b[0],
        a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
        a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0],
    ]


def carbon_isotopes(count: int) -> List[float]:
    """
    Calculates carbon isotope distributions.

    Args:
        count: Number of carbon atoms

    Returns:
        List of 4 floating point numbers representing isotope distributions
    """
    lambda_val = float(count) * 0.011
    c13 = [0.0] * 4
    fact = [1, 1, 2, 6]

    for k in range(4):
        c13[k] = pow(lambda_val, k) * np.exp(-lambda_val) / float(fact[k])

    return c13


def sulfur_isotopes(count: int) -> List[float]:
    """
    Calculates sulfur isotope distributions.

    Args:
        count: Number of sulfur atoms

    Returns:
        List of 4 floating point numbers representing convolved isotope distributions
    """
    lambda33 = float(count) * 0.0076
    lambda35 = float(count) * 0.044
    s33 = [0.0] * 4
    s35 = [
        pow(lambda35, 0) * np.exp(-lambda35),
        0.0,
        pow(lambda35, 1) * np.exp(-lambda35),
        0.0,
    ]

    fact = [1, 1, 2, 6]
    for k in range(4):
        s33[k] = pow(lambda33, k) * np.exp(-lambda33) / float(fact[k])

    return convolve(s33, s35)


def peptide_isotopes(carbons: int, sulfurs: int) -> List[float]:
    """
    Calculates peptide isotope distributions based on number of carbon and sulfur atoms.

    Args:
        carbons: Number of carbon atoms
        sulfurs: Number of sulfur atoms

    Returns:
        List of 3 floating point numbers representing normalized isotope distributions
    """
    c = carbon_isotopes(carbons)
    s = sulfur_isotopes(sulfurs)
    result = convolve(c, s)
    max_val = max(result[:3]).item()  # Only consider first 3 values for normalization

    # Normalize first 3 values
    return [val.item() / max_val for val in result[:3]]


def test_peptide_isotopes():
    """
    Test function to verify peptide isotope calculations.
    """
    iso = peptide_isotopes(60, 5)
    expected = [0.3972, 0.2824, 0.1869, 0.0846]
    expected = [val / 0.3972 for val in expected[:3]]  # Normalize first 3 values

    # Check if all differences are within tolerance
    tolerance = 0.02
    matched = all(abs(a - b) <= tolerance for a, b in zip(iso, expected))

    assert matched, f"Test failed: {iso} != {expected}"
    print("Test passed successfully!")


# Sample entry:

"""
 {
    "precursor": {
        "sequence": "PEPTIDEPINK",
        "charge": 2,
        "decoy": false
    },
    "elution_group": {
        "id": 0,
        "precursor_mzs": [
            1810.917339999999,
            1810.917339999999
        ],
        "fragment_mzs": {
            "a1": 123.0,
            "b1": 123.0,
            "c1^2": 123.0
        },
        "precursor_charge": 2,
        "mobility": 0.8,
        "rt_seconds": 0.0,
        "decoy": false,
        "expected_precursor_intensity": [
            1.0,
            1.0
        ],
        "expected_fragment_intensity": {
            "a1": 1.0,
            "b1": 1.0,
            "c1^2": 1.0
        }
    }
}
"""


def get_human_peptides():
    print("Cleaving the proteins with trypsin...")
    unique_peptides = set()
    fasta_file = "/Users/sebastianpaez/fasta/20231030_UP000005640_9606.fasta"
    # fasta_file = "/Users/sebastianpaez/git/timsseek/data/HeLa_cannonical_proteins.fasta"
    with open(fasta_file) as file:
        for description, sequence in fasta.FASTA(file):
            new_peptides = parser.cleave(
                sequence,
                "trypsin",
                min_length=6,
                max_length=20,
                missed_cleavages=0,
            )
            unique_peptides.update(new_peptides)

    unique_peptides = list(unique_peptides)
    print(unique_peptides[:5])
    print("Done, {0} sequences obtained!".format(len(unique_peptides)))
    return unique_peptides


def peptide_formula_dist(formula):
    c_count = 0
    s_count = 0

    for elem, _, count in formula.elements():
        str_elem = str(elem)
        if str_elem == "C":
            c_count = count
        elif str_elem == "S":
            s_count = count
        else:
            continue

    if c_count == 0:
        raise ValueError("No carbons found in formula")

    return peptide_isotopes(c_count, s_count)


def as_decoy(x: str) -> str:
    out = x[0] + x[-2:0:-1] + x[-1]
    return out


def yield_as_decoys(peptides: List[str]) -> Generator[str, None, None]:
    for peptide in peptides:
        yield as_decoy(peptide)


def yield_with_charges(
    peptides, min_charge, max_charge
) -> Generator[Tuple[str, int], None, None]:
    for peptide in peptides:
        for charge in range(min_charge, max_charge + 1):
            yield (peptide, charge)


PROTON_MASS = 1.007276466
NEUTRON_MASS = 1.008664916


def supersimpleprediction(mz, charge):
    intercept_ = -1.660e00
    log1p_mz = np.log1p(mz)
    sq_mz_over_charge = (mz**2) / charge
    log1p_sq_mz_over_charge = np.log1p(sq_mz_over_charge)

    out = (
        intercept_
        + (-3.798e-01 * log1p_mz)
        + (-2.389e-04 * mz)
        + (3.957e-01 * log1p_sq_mz_over_charge)
        + (4.157e-07 * sq_mz_over_charge)
        + (1.417e-01 * charge)
    )
    return out


def as_entry(
    peptide: rustyms.LinearPeptide,
    ion_dict: dict[str, (float, float)],
    decoy: bool,
    id: int,
) -> dict | None:
    pep_formula = peptide.formula()[0]
    isotope_dist = peptide_formula_dist(pep_formula)
    precursor_mzs = (
        pep_formula.mass() + (PROTON_MASS * peptide.charge)
    ) / peptide.charge
    if precursor_mzs < 400 or precursor_mzs > 1000:
        return None

    ims = float(supersimpleprediction(precursor_mzs, peptide.charge))

    neutron_fraction = NEUTRON_MASS / peptide.charge
    precursor_mzs = [
        float(precursor_mzs + (neutron_fraction * isotope)) for isotope in [-1, 0, 1, 2]
    ]

    max_intensity = max([v[1] for v in ion_dict.values()])
    max_keep = max_intensity * 0.02
    ion_dict = {
        k: v
        for k, v in ion_dict.items()
        if (v[0] > 250) and (v[1] > max_keep) and (v[0] < 2000)
    }

    ion_mzs = {k: v[0] for k, v in ion_dict.items()}
    ion_intensities = {k: v[1] for k, v in ion_dict.items()}
    if len(ion_mzs) < 3:
        return None

    precursor = {
        "sequence": str(peptide),
        "charge": peptide.charge,
        "decoy": decoy,
    }
    elution_group = {
        "id": id,
        "precursor_mzs": precursor_mzs,
        "fragment_mzs": ion_mzs,
        "precursor_charge": peptide.charge,
        "mobility": ims,
        "rt_seconds": 0.0,
        "decoy": decoy,
        "expected_precursor_intensity": [0.001] + isotope_dist,
        "expected_fragment_intensity": ion_intensities,
    }
    return {"precursor": precursor, "elution_group": elution_group}


def foo():
    # outfile = "FUUUUU_small.ndjson"
    outfile = "FUUUUU.ndjson"
    pretty_outfile = f"{outfile}.pretty.json"
    peps = get_human_peptides()
    decoys = list(yield_as_decoys(peps))

    pretty_outs = []
    is_first_n = 10
    id = 0

    pprint(f"Writing output to file: { outfile }")
    with open(outfile, "w") as file:
        targ_use = list(yield_with_charges(peps, 2, 3))
        for x in tqdm(
            model.predict_batched_annotated(
                targ_use,
                min_intensity=0.001,
                min_ordinal=3,
                max_ordinal=1000,
            ),
            desc="Targets",
            total=len(targ_use),
        ):
            id += 1
            elem = as_entry(x[0], x[1], False, id)
            if elem is None:
                continue
            if is_first_n > 0:
                pprint(elem)
                pretty_outs.append(elem)
                is_first_n -= 1

            file.write(json.dumps(elem) + "\n")
            file.flush()

        is_first_n = 10
        dec_use = list(yield_with_charges(decoys, 2, 3))
        for x in tqdm(
            model.predict_batched_annotated(
                dec_use,
                min_intensity=0.001,
                min_ordinal=3,
                max_ordinal=1000,
            ),
            desc="Decoys",
            total=len(dec_use),
        ):
            id += 1
            elem = as_entry(x[0], x[1], True, id)
            if elem is None:
                continue
            if is_first_n > 0:
                pprint(elem)
                pretty_outs.append(elem)
                is_first_n -= 1

            file.write(json.dumps(elem) + "\n")
            file.flush()

    pprint(f"Writing pretty output to file: { pretty_outfile }")
    with open(pretty_outfile, "w") as file:
        file.write(json.dumps(pretty_outs, indent=4))
        file.flush()


if __name__ == "__main__":
    # Run the test
    test_peptide_isotopes()
    foo()
