# Length-reweighted recovery gap bound

Published CATH-test-split recovery: 52.4%. Medium-length (70-150 aa) structures are the closest available proxy to typical CATH-domain sizes; comparing recovery on that subset alone to the full 30-structure mean bounds how much of the gap is length-distribution skew rather than genuine train/test mismatch.

| model | all-30 mean | medium-only mean | shift toward published | gap all-30 | gap medium-only |
|---|---|---|---|---|---|
| fullbackbone | 46.4% | 49.0% | +2.6 pp | +6.0 pp | +3.4 pp |
| ca_only | 41.3% | 41.9% | +0.6 pp | +11.1 pp | +10.5 pp |

(full-backbone is the model paper's headline comparison; published number is full-backbone T=0.1.)
