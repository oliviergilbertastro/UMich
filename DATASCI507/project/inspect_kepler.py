from collections import Counter
from datasets import load_dataset

tces = load_dataset(
    "nyu-planet-finders/kepler-raw-lightcurves",
    "tces",
    split="full",
    streaming=True,
)

labels = Counter()
star_ids = set()
rows = 0

for tce in tces:
    labels[tce["av_training_set"]] += 1
    star_ids.add(tce["kepid"])
    rows += 1

print("TCE rows:", rows)
print("Unique stars:", len(star_ids))
print("Rows by label:", labels)