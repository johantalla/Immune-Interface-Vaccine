import numpy as np
import matplotlib.pyplot as plt

x = np.load(
    "./boltz_outputs/boltz_results_boltz_inputs/predictions/wzg_seq29/"
    "plddt_wzg_seq29_model_0.npz"
)

plddt = x["plddt"] * 100

plt.figure(figsize=(12, 4))

plt.plot(
    range(1, len(plddt) + 1),
    plddt
)

plt.xlabel("Residue position")
plt.ylabel("pLDDT")
plt.ylim(0, 100)

plt.title("Boltz-2 prediction confidence for Wzg seq29")

plt.tight_layout()
plt.show()