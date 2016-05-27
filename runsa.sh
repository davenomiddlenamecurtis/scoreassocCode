if [ .$1 == . ]
then
echo Usage: runsa.sh GENE
echo
echo This is a simple script to demonstrate running scoreassoc using plink-seq files
echo First set up a plink-seq project (see plink-seq documentation)
echo Then edit the script to set your project name and your phenotype
echo You need to have files called funcWeights.txt and exclusions.txt in the current working directory
echo (or else edit the script to point to their correct locations)
exit 
fi
		
PROJECT=SSS
PHENO=scz
# replace above lines with your project and phenotype names

GENE=$1

pseq ${PROJECT} v-view --geno  --phenotype ${PHENO} --gene $GENE > ${GENE}.dat
pseq ${PROJECT} counts --annotate refseq --gene $GENE --phenotype ${PHENO} > ${GENE}.annot
pscoreassoc ${GENE}.dat --annotfile ${GENE}.annot --weightfile funcWeights.txt --filterfile exclusions.txt --outfile ${GENE}.psao --minweight 80 --ldthreshold 0.7 --dorecessive

