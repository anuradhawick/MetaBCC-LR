import math
import numpy as np


def _line_to_vec(line):
    return np.fromstring(line.strip(), sep=" ", dtype=float)


def _load_bins(stats_file):
    with open(stats_file, "r") as handle:
        lines = [line.rstrip("\n") for line in handle]

    bins = []
    i = 0
    while i < len(lines):
        name = lines[i].strip()
        i += 1
        if name == "":
            continue
        if i + 3 >= len(lines):
            raise ValueError("Invalid cluster-stats format: incomplete bin statistics block.")

        p15_mean = _line_to_vec(lines[i])
        p3_mean = _line_to_vec(lines[i + 1])
        p15_std = _line_to_vec(lines[i + 2])
        p3_std = _line_to_vec(lines[i + 3])
        i += 4

        bins.append(
            {
                "name": name,
                "p15_mean": p15_mean,
                "p15_std": p15_std,
                "p3_mean": p3_mean,
                "p3_std": p3_std,
            }
        )

    if len(bins) == 0:
        raise ValueError("No bins found in cluster-stats file.")

    return bins


def _normal_log_prob(x, mu, std):
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (x - mu) / std
        return (-0.5 * (z ** 2) - np.log(math.sqrt(2.0 * math.pi) * std)).sum(axis=1)


def assign_reads_numpy(p3_path, p15_path, stats_file, output_path, batch_size=100000):
    bins = _load_bins(stats_file)

    with open(p3_path, "r") as p3_handle, open(p15_path, "r") as p15_handle, open(output_path, "w+") as output:
        batch3 = []
        batch15 = []

        for p3_line, p15_line in zip(p3_handle, p15_handle):
            batch3.append(_line_to_vec(p3_line))
            batch15.append(_line_to_vec(p15_line))

            if len(batch3) >= batch_size:
                _write_batch(batch3, batch15, bins, output)
                batch3.clear()
                batch15.clear()

        if len(batch3) > 0:
            _write_batch(batch3, batch15, bins, output)


def _write_batch(batch3, batch15, bins, output_handle):
    x3 = np.vstack(batch3)
    x15 = np.vstack(batch15)

    n = x3.shape[0]
    best_prob15 = np.full(n, -np.inf, dtype=float)
    best_prob3 = np.full(n, -np.inf, dtype=float)
    best_names = np.full(n, "UnBinned", dtype=object)

    for b in bins:
        prob15 = _normal_log_prob(x15, b["p15_mean"], b["p15_std"])
        prob3 = _normal_log_prob(x3, b["p3_mean"], b["p3_std"])
        better = (prob15 > best_prob15) | ((prob15 == best_prob15) & (prob3 > best_prob3))
        best_prob15[better] = prob15[better]
        best_prob3[better] = prob3[better]
        best_names[better] = b["name"]

    output_handle.write("\n".join(best_names.tolist()) + "\n")
