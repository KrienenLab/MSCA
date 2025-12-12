#!/usr/bin/env bash
# collect summary.csv files into a single csv file
cd $(dirname $0)
touch cellranger-arc_summaries.csv
cat sample_names.txt | while read SAMPLE; do
    SAMPLE=$(echo $SAMPLE | sed 's:/$::')
    [ -e $SAMPLE/$SAMPLE/outs/summary.csv ] || continue
    HEADER=$(head -1 "$SAMPLE/$SAMPLE/outs/summary.csv")
    DATA=$(grep -v "$HEADER" "$SAMPLE/$SAMPLE/outs/summary.csv")
    # if file is empty, add the header
    if [ ! -s cellranger-arc_summaries.csv ]; then
        echo "$HEADER" >> cellranger-arc_summaries.csv
    fi
    echo "$DATA" >> cellranger-arc_summaries.csv
done
exit 0
