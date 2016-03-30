#include "scoreassoc.hpp"

#include "dcexpr.hpp"
#include "safilterfuncs.hpp"

#define MAXEXCLUSIONS 20
#define MAXEXCLUSIONSTRLENGTH 200

int currentLocus,nExc;

char exclusionStr[MAXEXCLUSIONS][MAXEXCLUSIONSTRLENGTH];

express exclusion[MAXEXCLUSIONS];
#define EVAL_R1 \
if ((r1=b1->eval())==NULL) return NULL; 
#define EVAL_R1_R2 \
if ((r1=b1->eval())==NULL || (r2=b2->eval())==NULL) return NULL; 

dcexpr_val *alt_count_op(vnode *b1)
{
dcexpr_val *r1;
int s,i,c;
EVAL_R1;
s=double(*r1);
double rv=0;
for (i=0;i<2;++i)
//	if (global_sub[s-1]->all[currentLocus][i]==rarer[currentLocus]) rarer not set yet
	if (global_sub[s-1]->all[currentLocus][i]==2)
		++rv;
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *weight_op(vnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=weight[currentLocus];
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *freq_op(vnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=cc_freq[f][currentLocus];
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *geno_count_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
int cc,g;
EVAL_R1_R2;
cc=double(*r1);
g=double(*r2);
double rv=cc_genocount[cc][g][currentLocus];
delete r1;
delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *nsub_op(vnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=cc_count[f][currentLocus];
delete r1;
return new dcexpr_double(rv);
}

int stateExclusions(FILE *fp)
{
	int e;
	if (nExc!=0)
	{
		fprintf(fp,"\nNB Any loci meeting the following exclusion criteri%s were not included in the above analysis:\n",nExc==1?"on":"a");
		for (e=0;e<nExc;++e)
			fprintf(fp,"%s\n",exclusionStr[e]);
	}
	return 1;
}

int initExclusions(FILE *fp,char *extras[])
{
	int e;
	add_un_op("FREQ",freq_op); // frequency of this variant in controls or cases (0 or 1)
	add_un_op("NSUB",nsub_op); // number of control or case subjects (0 or 1) typed for this variant 
	add_un_op("ALTCOUNT",alt_count_op);	// number of ALT alleles possessed by subject s
	add_un_op("WEIGHT",weight_op);	// weight of variant - this is the supplied functional weight because exclusions applied before frequency weights are calculated
	add_bin_op_next("GENOCOUNT",geno_count_op); // number of case or control subjects who are have genototype - cc GENOCOUNT g AA (0 or 1, 0, 1 or 2)
	for (nExc=0;fgets(long_line,LONG_LINE_LENGTH,fp) && sscanf(long_line,"%[^\n]",exclusionStr[nExc])==1;++nExc)
	{
		assert(nExc<MAXEXCLUSIONS);
		if (exclusion[nExc].parse(exclusionStr[nExc])==0)
			return 0;
	}
	if (extras!=0)
		for (e=0;extras[e][0]!='\0';++e,nExc++)
		{
			strcpy(exclusionStr[nExc],extras[e]);
			assert(nExc<MAXEXCLUSIONS);
			if (exclusion[nExc].parse(exclusionStr[nExc])==0)
				return 0;
		}
	return 1;
}

int applyExclusions(par_info *pi)
{
	int l,ll,e;
for (l=0;l<pi->n_loci_to_use;++l)
{
	currentLocus=pi->loci_to_use[l];
	for (e=0;e<nExc;++e)
	{
		if (double(*exclusion[e].eval())!=0.0)
		{
			printf("Excluding locus %d (%s) because it fails this condition: %s(contfreq=%f,casefreq=%f)\n",currentLocus+1,names[currentLocus],exclusionStr[e],cc_freq[0][currentLocus],cc_freq[1][currentLocus]);
			for (ll=l;ll<pi->n_loci_to_use-1;++ll)
				pi->loci_to_use[ll]=pi->loci_to_use[ll+1];
			--pi->n_loci_to_use;
			--l;
			break;
		}
	}
}
	return 1;
}
