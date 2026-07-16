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

awk 'NR==FNR && NR>1 {pop[$1]=$2; next}
NR>1 {
  kin=$8
  if      (kin > 0.354)          rel="Duplicate_or_MZ";
  else if (kin >= 0.177)         rel="1st";
  else if (kin >= 0.0884)        rel="2nd";
  else if (kin >= 0.0442)        rel="3rd";
  else                           rel="Unrel";

  id1=$2; id2=$4;

  deg[id1,rel]=deg[id1,rel] ? deg[id1,rel]","id2 : id2;
  deg[id2,rel]=deg[id2,rel] ? deg[id2,rel]","id1 : id1;
}
END {
  print "IID\tpop\thas_1st_deg\thas_2nd_deg\thas_3rd_deg\tfirst_ids\tsecond_ids\tthird_ids";

  for (id in pop) {
    first = (deg[id,"1st"] != "") ? "Yes" : "No";
    second = (deg[id,"2nd"] != "") ? "Yes" : "No";
    third = (deg[id,"3rd"] != "") ? "Yes" : "No";

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
           id, pop[id], first, second, third,
           (deg[id,"1st"]  ? deg[id,"1st"]  : "-"),
           (deg[id,"2nd"]  ? deg[id,"2nd"]  : "-"),
           (deg[id,"3rd"]  ? deg[id,"3rd"]  : "-");
  }
}' depth.tsv Diacarnus0526_king.kin0 > king_individual_summary.tsv
