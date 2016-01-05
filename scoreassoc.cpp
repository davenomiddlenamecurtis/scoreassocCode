/* scoreassoc.cpp */
extern "C"
{
#include "scoreassoc.h"
};

#include <ctype.h>

#ifndef NOFILTERS
#include "dcexpr.hpp"
#include "safilterfuncs.hpp"
#endif

#define PROGRAM "SCOREASSOC"
#define SAVERSION "3.1"

#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

int main(int argc,char *argv[])
{
FILE *fp,*fd,*fo;
int nsub,l,n_new_sub,real_nsub;
float *score;
int max_cc[2],s,n_non_mendelian;
non_mendelian *non_mendelians;
char *non_mendelian_report;
par_info pi;
sa_par_info spi;
subject **sub,**new_sub,**real_sub;
float p;
printf("%s v%s\n",PROGRAM,SAVERSION);
printf("MAX_LOCI=%d\nMAX_ALL=%d\nMAX_SUB=%d\n",MAX_LOCI,
	   MAX_ALL,MAX_SUB);

assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
for (s=0;s<MAX_SUB;++s)
	assert(sub[s]=(subject *)calloc(1,sizeof(subject)));
assert(score=(float *)calloc(MAX_SUB,sizeof(float)));
assert(argc>=4);
assert((fp=fopen(argv[1],"r"))!=0);
assert((fd=fopen(argv[2],"r"))!=0);
assert((fo=fopen(argv[3],"w"))!=0);
fprintf(fo,"SCOREASSOC output\n");
read_par(fp,&pi);
for (l=0;l<pi.n_loci_to_use;++l)
  if (pi.n_alls[pi.loci_to_use[l]]!=2)
    { error("All loci used must be biallelic",""); return 1; }
read_sa_par(fp,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
//read_score_assoc_par(fp,&pi,&wfactor,&use_func_weights,use_cc_freqs,func_weight,cc_freq,cc_count,max_cc,&use_locus_names,names,&use_comments,comments,&do_recessive_test,&weight_threshold,&LD_threshold,&use_haplotypes);
if (pi.use_cc==0 || pi.use_cc==2)
  { error("Need to use case-control data",""); return 1; }
if (spi.use_trios)
{
	if (atoi(comments[0])>22 || toupper(comments[0][0]) == 'X' || toupper(comments[0][0]) == 'Y' ||
		toupper(comments[0][0]) == 'C'&&toupper(comments[0][1]) == 'H'&&toupper(comments[0][2]) == 'R' && (atoi(comments[0] + 3) > 22 || toupper(comments[0][3]) == 'X' || toupper(comments[0][3]) == 'Y'))
	{
		error("Cannot at present use trios for genes on X or Y chromosome", "");
		return 1;
	}
}

global_sub=sub;
initExclusions(fp);
fclose(fp);
if (spi.use_cc_freqs[0]==0 || spi.use_cc_freqs[1]==0)
	read_all_subjects(fd,sub,&nsub,&pi);
else
	nsub=0;
fclose(fd);
if (spi.use_trios)
	{
		assert(new_sub=(subject **)calloc(nsub,sizeof(subject*)));
		for (s=0;s<nsub;++s)
			assert(new_sub[s]=(subject *)calloc(1,sizeof(subject)));
		assert(non_mendelians=(non_mendelian *)calloc(MAX_SUB,sizeof(non_mendelian)));
		if ((n_new_sub=sort_trios(sub,nsub,&pi,new_sub,non_mendelians,&n_non_mendelian,long_line))==0)
			exit(1);
		real_sub=sub;
		sub=new_sub;
		real_nsub=nsub;
		nsub=n_new_sub;
		assert(non_mendelian_report=(char*)malloc(strlen(long_line)+1));
		strcpy(non_mendelian_report,long_line);
	}
fprintf(fo,
"Locus                   controls     frequency        cases          frequency   frequency allele  weight\n"
"                     AA  :   AB  :  BB                  AA  :   AB  :  BB                     \n");
get_freqs(sub,nsub,&pi,&spi,cc_freq,cc_count,cc_genocount);
applyExclusions(&pi);
set_weights(fo,weight,missing_score,rarer,sub,nsub,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
get_scores(score,weight,missing_score,rarer,sub,nsub,&pi);
if (argc>4)
	write_scores(argv[4],sub,nsub,score);

p=do_score_onetailed_ttest(fo,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer);
#ifdef DO_WILCOXON
if (spi.use_cc_freqs[0]==0 && spi.use_cc_freqs[1]==0)
	do_score_onetailed_wilcoxon(fo,score,sub,nsub);
#endif
// p=do_score_recessive_test(fo,score,sub,nsub,&pi,use_cc_freqs,cc_freq,cc_count,max_cc,weight,missing_score,rarer);
if (spi.do_recessive_test)
{
if (atoi(comments[0])>22 || toupper(comments[0][0])=='X' || toupper(comments[0][0])=='Y' || 
	toupper(comments[0][0])=='C'&&toupper(comments[0][1])=='H'&&toupper(comments[0][2])=='R'&&(atoi(comments[0]+3)>22 || toupper(comments[0][3])=='X' || toupper(comments[0][3])=='Y'))
// simple trick for now to avoid X and Y genes
	fprintf(fo,"\nCannot do recessive test for genes on X or Y chromosome.\n");
else if (spi.use_cc_freqs[0] || spi.use_cc_freqs[1])
	fprintf(fo,"\nCannot do recessive test with supplied allele frequencies instead of genotypes\n");
else if (spi.use_haplotypes)
	do_recessive_HWE_test_with_haplotypes(fo,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
else 
	do_recessive_HWE_test(fo,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
}

if (spi.use_trios)
{
	fprintf(fo,"\n%s\n",*non_mendelian_report?non_mendelian_report:"No non_mendelian transmissions or de novo mutations found\n");
}

#ifdef USEFILTERS
stateExclusions(fo);
#endif
fclose(fo);
printf("\nProgram run completed OK\n");
return 0;
}
 