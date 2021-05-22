import random
import numpy as np
import logging

logger = logging.getLogger('MetaBCC-LR')

def sample(output, sample_count, ground_truth):
    p3_data = np.loadtxt(f"./{output}/profiles/3mers", dtype=float)
    p15_data = np.loadtxt(f"./{output}/profiles/15mers", dtype=float)

    logger.debug(f"3mer data shape {str(p3_data.shape)}")
    logger.debug(f"15mer data shape {str(p15_data.shape)}")

    if ground_truth is not None:
        ground_truth = np.array(open(ground_truth).read().strip().split('\n'))
        logger.debug(f"Ground truth data shape {str(ground_truth.shape)}")

    if sample_count <= 0:
        sample_count = int(len(p3_data)/100)
    
    logger.debug(f"Sampling count {sample_count}")

    idx_chosen = random.sample(range(len(p3_data)), sample_count)

    p3_sampled = p3_data[idx_chosen]
    p15_sampled = p15_data[idx_chosen]

    np.save(f"./{output}/profiles/3mers_sampled.npy", p3_sampled)
    np.save(f"./{output}/profiles/15mers_sampled.npy", p15_sampled)

    if ground_truth is not None:
        ground_truth_sampled = ground_truth[idx_chosen]
        np.save(f"./{output}/misc/filtered_truth_sampled.npy", ground_truth_sampled)