for method in avg_clade leaf_dist_avg leaf_dist_max leaf_dist_min length length_clade max max_clade med_clade root_dist single_linkage single_linkage_cut single_linkage_union sum_branch sum_branch_clade
do
    cat thresholds.txt | parallel "TreeCluster.py -t {} -m $method -o {}.${method}.t.txt -i No_clade_13_phylogeny.txt"
done

