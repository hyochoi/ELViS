cat package_conflict_warnning.txt | grep -Po '[a-zA-Z.]+::[a-zA-Z._]+' | sed 's/::/\t/g' | sort | uniq > package_conflict_warnning.txt_parsed.txts

grep -Po $(cut -f 2 package_conflict_warnning.txt_parsed.txt | sort | uniq | awk  'NR!=1{printf"|"}{printf "[^[:space:]]+"$1}') ../R/*  | sort |
uniq > package_conflict_warnning.txt_parsed.txt_to_review.txt


echo \
'base::intersect
dplyr::slice
stats::median
stats::quantile
generics::setdiff' | sed 's/::/\t/g' | sort | uniq > package_conflict_warnning.txt_parsed.txt_to_review.txt_fns_to_keep.txt

cat package_conflict_warnning.txt_parsed.txt package_conflict_warnning.txt_parsed.txt_to_review.txt_fns_to_keep.txt  | sort | uniq | cat /dev/stdin package_conflict_warnning.txt_parsed.txt_to_review.txt_fns_to_keep.txt | sort | uniq -c | awk '$1==1' > package_conflict_warnning.txt_parsed.txt_fns_to_exclude.tx


cat package_conflict_warnning.txt_parsed.txt_fns_to_exclude.txt |
awk '{a[$2]=a[$2]","$3}END{for(i in a) { sub("^,","",a[i]);print "#'"'"' @rawNamespace import("i", except=c("a[i]"))" }}' \
 > package_conflict_warnning.txt_parsed.txt_fns_to_exclude.txt_roxygen2_import_strings.txt


grep -P $(cat package_conflict_warnning.txt_parsed.txt_fns_to_exclude.txt | awk '{print $2}' | sort | uniq | awk 'NR!=1{printf "|"}{printf $1}' ) ../R/* | grep @import \
 > package_conflict_warnning.txt_parsed.txt_fns_to_exclude.txt_which_files.txt
