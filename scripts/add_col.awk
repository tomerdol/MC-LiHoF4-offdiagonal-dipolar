# script used to add new columns "mk2y","mk2z","tMaxConfigs" to existing results files
# used by add_col_to_results.sh
{ if ($1 ~ /^#/ && $2 != "index") {print $0}
else if ($1 ~ /^#/ && $2 == "index") {printf "#%10s %17s %19s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %5s\n", $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,"mk2y","mk2z","tMaxConfigs",$17}
else if ($1 == "index") {printf " %10s %17s %19s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %5s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"mk2y","mk2z","tMaxConfigs",$16}
else {printf " %10s %17s %19s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %17s %5s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"NaN","NaN","NaN",$16}}
