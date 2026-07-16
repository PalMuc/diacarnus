cat > depth.tsv << 'EOF'
IID	pop
GW34916	Shallow
GW34917	Shallow
GW34918	Shallow
GW34919	Shallow
GW34920	Shallow
GW34921	Shallow
GW34922	Shallow
GW34923	Shallow
GW34924	Shallow
GW34925	Shallow
GW34926	Shallow
GW34927	Shallow
GW34928	Mesophotic
GW34929	Mesophotic
GW34930	Mesophotic
GW34931	Mesophotic
GW34932	Mesophotic
GW34933	Mesophotic
GW34934	Mesophotic
EOF

awk 'NR==FNR && NR>1 {pop[$1]=$2; next}
NR>1 {
  kin=$8
  if      (kin > 0.354)          rel="Duplicate_or_MZ";
  else if (kin >= 0.177)         rel="1st_degree";
  else if (kin >= 0.0884)        rel="2nd_degree";
  else if (kin >= 0.0442)        rel="3rd_degree";
  else                           rel="Unrelated";

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
         $2, pop[$2], $4, pop[$4], kin, rel, $6, $7
}' depth.tsv Diacarnus0526_king.kin0 > king_pairs_annotated.tsv

###Then on python###

#!/usr/bin/env python3
import csv

pairs_file = "king_pairs_annotated.tsv"
depth_file = "depth.tsv"
out_file = "king_individual_summary.tsv"

# Load depth info, using first column as IID and second as pop
id2pop = {}
with open(depth_file) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)  # assume there is a header
    # first column = IID, second column = pop (whatever its name is)
    for row in reader:
        if not row or len(row) < 2:
            continue
        iid = row[0].strip()
        pop = row[1].strip()
        id2pop[iid] = pop

# Initialize per-individual containers
rel_per_ind = {
    iid: {
        "1st_degree": set(),
        "2nd_degree": set(),
        "3rd_degree": set()
    }
    for iid in id2pop
}

# Read pairwise relationships from king_pairs_annotated.tsv
with open(pairs_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        iid1 = row["IID1"].strip()
        iid2 = row["IID2"].strip()
        rel = row["Duplicate_or_MZ"].strip()

        if rel in {"1st_degree", "2nd_degree", "3rd_degree"}:
            rel_per_ind[iid1][rel].add(iid2)
            rel_per_ind[iid2][rel].add(iid1)

# Write per-individual summary
with open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(
        ["IID", "pop", "has_1st_deg", "has_2nd_deg", "has_3rd_deg",
         "first_ids", "second_ids", "third_ids"]
    )
    for iid, pop in id2pop.items():
        first = sorted(rel_per_ind[iid]["1st_degree"])
        second = sorted(rel_per_ind[iid]["2nd_degree"])
        third = sorted(rel_per_ind[iid]["3rd_degree"])
        writer.writerow([
            iid,
            pop,
            "Yes" if first else "No",
            "Yes" if second else "No",
            "Yes" if third else "No",
            ",".join(first) if first else "-",
            ",".join(second) if second else "-",
            ",".join(third) if third else "-"
        ])
