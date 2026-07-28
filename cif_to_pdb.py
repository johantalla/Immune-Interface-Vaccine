import sys
import os
import glob
import logging
from Bio.PDB import PDBIO, MMCIFParser

class OutOfChainsError(Exception):
    pass

def int_to_chain(i, base=62):
    """Convert an integer to a chain ID (A-Z, a-z, 0-9)"""
    chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if i < base:
        return chars[i]
    raise OutOfChainsError()

def rename_chains(structure):
    """Renames chains to be one-letter chains, returns old->new mapping"""
    next_chain = 0
    chainmap = {c.id: c.id for c in structure.get_chains() if len(c.id) == 1}
    for o in structure.get_chains():
        if len(o.id) != 1:
            if o.id[0] not in chainmap:
                chainmap[o.id[0]] = o.id
                o.id = o.id[0]
            else:
                c = int_to_chain(next_chain)
                while c in chainmap:
                    next_chain += 1
                    c = int_to_chain(next_chain)
                    if next_chain >= 62:
                        raise OutOfChainsError()
                chainmap[c] = o.id
                o.id = c
    return chainmap

def convert_cif_to_pdb(ciffile):
    pdbfile = ciffile.replace(".cif", ".pdb")
    parser = MMCIFParser(QUIET=True)
    strucid = os.path.basename(ciffile)[:4] if len(os.path.basename(ciffile)) > 4 else "1xxx"
    structure = parser.get_structure(strucid, ciffile)

    try:
        rename_chains(structure)
    except OutOfChainsError:
        logging.error(f"Too many chains to represent in PDB format: {ciffile}")
        return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdbfile)
    return True

# Batch process every .cif file
cif_files = glob.glob("boltz_outputs/boltz_results_boltz_inputs/predictions/**/*.cif", recursive=True)
os.makedirs("pdb_outputs", exist_ok=True)

failed = []
for i, ciffile in enumerate(cif_files):
    seq_name = os.path.basename(os.path.dirname(ciffile))
    out_path = os.path.join("pdb_outputs", f"{seq_name}.pdb")
    try:
        success = convert_cif_to_pdb(ciffile)
        if success:
            # move/rename to pdb_outputs with the sequence name
            generated_pdb = ciffile.replace(".cif", ".pdb")
            os.rename(generated_pdb, out_path)
        else:
            failed.append(ciffile)
    except Exception as e:
        failed.append((ciffile, str(e)))
    if (i + 1) % 100 == 0:
        print(f"Converted {i + 1}/{len(cif_files)}")

print(f"\nDone. {len(cif_files) - len(failed)} converted, {len(failed)} failed.")
if failed:
    print("Failed files:", failed)