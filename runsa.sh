if [ .$1 == . ]
then
echo Usage: runpscore.sh GENE
exit
fi
		
PROJECT=SSS
PHENO=scz
GENE=$1
pseq ${PROJECT} v-view --geno  --phenotype ${PHENO} --gene $GENE > ${GENE}.dat
pseq ${PROJECT} counts --annotate refseq --gene $GENE --phenotype ${PHENO} > ${GENE}.annot
pscoreassoc ${GENE}.dat --annotfile ${GENE}.annot --weightfile funcWeights.txt --filterfile exclusions.txt --outfile ${GENE}.psao --minweight 80 --ldthreshold 0.7 --dorecessive

