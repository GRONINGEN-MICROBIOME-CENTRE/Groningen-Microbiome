#{input.Polymor} [input.Mutation} {input.mutations_c} {ouput.poly_c} {output.cov} {input.Covered_positions}
Mutations_I=$2
Mutations_P=$1

OUT1=$3
OUT2=$4
OUT3=$5

FILE=$6



awk 'BEGIN{FS=OFS="\t"}NR==FNR{z[$1"|"$2]=$2;next}{if (z[$1"|"$2]){print $0"\tMutation"}}' $FILE $Mutations_I > $OUT1

awk 'BEGIN{FS=OFS="\t"}NR==FNR{z[$1"|"$2"]=$2;next}
{if (z[$1"|"$2]){print $0"\tPolymorphism"}}' $FILE $Mutations_P  > "$OUT2"
cat $OUT1 $OUT2 >  "$OUT3"
