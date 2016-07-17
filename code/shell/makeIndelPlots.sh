if [[ "$1" == "" || "$2" == "" || "$3" == "" || "$4" == "" || "$5" == "" ]] ; then
    printf "Usage: source code/shell/makeIndelPlots.sh <baseFileName> <near> <maxIndelLen> <featureRE> <speciesPair>\n"
    printf "Example: source code/shell/makeIndelPlots.sh outCL13/IndelGroupsNonoverlappingWithInNearFeatures_K13k4L100D1_3000A100_3000d100_100N2F0 1000 3000 '^[^-]+-([^-(]+)[-(]' Col-0,Ler-0\n"
    printf "Example: source code/shell/makeIndelPlots.sh outHP14/IndelGroupsNonoverlappingWithInNearFeatures_K14k4L100D1_3000A100_3000d100_100N2F0 1000 3000 '(^[^:]+):' S.lycopersicum,S.pennellii\n"
    return 1
fi

Rscript code/R/plotIndelsWithFeatures.R . $1.indels.tsv $1.indels.pdf $2 $3 $4 $5
