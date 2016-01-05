#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string>
#include <map>

extern "C"
{
#include "scoreassoc.h"
};

#include "dcexpr.hpp"
#include "safilterfuncs.hpp"

#define PROGRAM "pscoreassoc"
#define PSAVERSION "3.1"

/*
Use same formats as pseq
 - same phenotype files
 - output of locus files with annotations and genotypes
 Use -- for all options
 Filters in their own text file
 Other options on command line
 Maintain compatibility with scoreassoc
 http://atgu.mgh.harvard.edu/plinkseq/dist/version-0.08/plinkseq-0.08-x86_64.tar.gz

 pseq SSS v-view --geno  --phenotype scz --gene $GENE > ${GENE}.dat
 pseq SSS counts --mask --annotate refseq --gene $GENE > ${GENE}.annot


 */

enum OPT {
	WEIGHTFACTOR=0,OUTFILE,WEIGHTFILE,ANNOTFILE,FILTERFILE,LDTHRESHOLD,WEIGHTTHRESHOLD,DORECESSIVE,USEHAPS,TRIOFILE,NUMOPTS
};

struct option_t {
	char *str;
	OPT o;
};

typedef struct option_t option;

option opt[]=
{
	{ "weightfactor", WEIGHTFACTOR },
	{ "outfile", OUTFILE },
	{ "weightfile", WEIGHTFILE },
	{ "annotfile", ANNOTFILE },
	{ "filterfile", FILTERFILE },
	{ "ldthreshold", LDTHRESHOLD },
	{ "minweight", WEIGHTTHRESHOLD },
	{ "dorecessive", DORECESSIVE },
	{ "usehaps", USEHAPS },
	{ "triofile", TRIOFILE },
	{"bad option", NUMOPTS}
};

class psa_par_info {
	// additional options for pscoreassoc not used by scoreassoc
public:
	psa_par_info();
	~psa_par_info();
	char outfile_fn[200],weightfile_fn[200],datafile_fn[200],filterfile_fn[200],annotfile_fn[200];
	FILE *outfile,*weightfile,*datafile,*filterfile,*annotfile;
};

psa_par_info::psa_par_info()
{
	outfile_fn[0]=weightfile_fn[0]=datafile_fn[0]=filterfile_fn[0]=annotfile_fn[0]='\0';
	outfile=stdout;
	weightfile=datafile=filterfile=annotfile=0;
}

psa_par_info::~psa_par_info()
{
	if (outfile!=stdout && outfile)
		fclose(outfile);
	if (weightfile)
		fclose(weightfile);
	if (datafile)
		fclose(datafile);
	if (filterfile)
		fclose(filterfile);
	if (annotfile)
		fclose(annotfile);
}


void usage()
{
	printf("pscoreassoc datafile [options]\n\nOptions:\n"
		"--weightfactor x (weight for very rare variants, default 10)\n"
"--outfile file\n"
"--weightfile file (specify functional weights)\n"
"--annotfile file (annotations from plink/seq)\n"
"--filterfile file\n"
"--ldthreshold x (to discard variants in LD for recessive analysis, default 0.9)\n"
"--minweight x (to include in recessive analysis, default 0, i.e. all variants)\n"
"--dorecessive\n"
"--usehaps\n"
"--triofile file\n"
);
}

char *skip_word(char *ptr)
{
	while (*ptr && !isspace(*ptr))
		++ptr;
	while (*ptr && isspace(*ptr))
		++ptr;
	return ptr;
}

int make_arg_string(char *arg_string, int argc, char *argv[])
{
	// one day might get some of these from a file rather than just from command line
	int a;
	*arg_string='\0';
	for (a = 1; a < argc; ++a)
	{
		strcat(arg_string,argv[a]);
		strcat(arg_string," ");
	}
	return 1;
}

int parse_arg_string(char *args, par_info *pi, sa_par_info *spi, psa_par_info *pspi)
{
	char *ptr;
	int a;
	spi->use_locus_names=spi->use_comments=1;
	spi->wfactor=10;
	spi->do_recessive_test=0;
	spi->LD_threshold=0.9;
	spi->use_haplotypes=spi->use_trios=0;
	spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=0;
	if (sscanf(args, "%s", pspi->datafile_fn) != 1 || pspi->datafile_fn[0]=='-')
	{
		usage();
		exit(1);
	}
	ptr=args;
	ptr=skip_word(ptr);
	while (*ptr)
	{
		if (strncmp(ptr, "--", 2))
		{
			printf("Need to specify options with --\n%s\n\n", ptr);
			usage();
			exit(1);
		}
		ptr+=2;
		for (a=0;a<NUMOPTS;++a)
			if (!strncmp(ptr,opt[a].str,strlen(opt[a].str)))
				break;
		if (a == NUMOPTS)
			{
			printf("Unrecognised option: --\n%s\n\n", ptr);
			usage();
			exit(1);
			}
		ptr=skip_word(ptr);
		int error=0,eatarg=0;
		switch (opt[a].o)
		{
		case WEIGHTFACTOR:
			if (sscanf(ptr, "%f", &spi->wfactor) != 1)
				error=1;
			else
				eatarg=1;
			break;
		case OUTFILE:
			if (*ptr=='-' || sscanf(ptr, "%s", pspi->outfile_fn) != 1)
				error=1;
			else
				eatarg=1;
			break;
		case WEIGHTFILE:
			if (*ptr=='-' || sscanf(ptr, "%s", pspi->weightfile_fn) != 1)
				error=1;
			else
				eatarg=1;
			break;
		case ANNOTFILE:
			if (*ptr=='-' || sscanf(ptr, "%s", pspi->annotfile_fn) != 1)
				error=1;
			else
				eatarg=1;
			break;
		case TRIOFILE:
			if (*ptr=='-' || sscanf(ptr, "%s", trios_fn) != 1)
				error=1;
			else
			{
				spi->use_trios=1;
				eatarg = 1;
			}
			break;
		case FILTERFILE:
			if (*ptr=='-' || sscanf(ptr, "%s", pspi->filterfile_fn) != 1)
				error=1;
			else
				eatarg=1;
			break;
		case LDTHRESHOLD:
			if (sscanf(ptr,"%f",&spi->LD_threshold)!=1)
				error=1;
			else
				eatarg=1;
			break;
		case WEIGHTTHRESHOLD:
			if (sscanf(ptr,"%f",&spi->weight_threshold)!=1)
				error=1;
			else
				eatarg=1;
			break;
		case DORECESSIVE:
			spi->do_recessive_test=1;
			break;
		case USEHAPS:
			spi->use_haplotypes=1;
			break;
		}
		if (error)
			{
				dcerror(1, "Error reading %s\n %s", opt[a].str,ptr); exit(1);
			}
		if (eatarg)
			ptr=skip_word(ptr);
		else 
			while (*ptr && isspace(*ptr))
				++ptr;
	}
	return 1;
}

void process_options(par_info *pi, sa_par_info *spi, psa_par_info *pspi)
{
	if (pspi->outfile_fn[0] && (pspi->outfile = fopen(pspi->outfile_fn, "w")) == 0)
	{
		dcerror(1, "Could not open outfile %s\n", pspi->outfile_fn); exit(1);
	}
	if (pspi->weightfile_fn[0] && (pspi->weightfile = fopen(pspi->weightfile_fn, "r")) == 0)
	{
		dcerror(1, "Could not open weightfile %s\n", pspi->weightfile_fn); exit(1);
	}
	if (pspi->annotfile_fn[0] && (pspi->annotfile = fopen(pspi->annotfile_fn, "r")) == 0)
	{
		dcerror(1, "Could not open annotfile %s\n", pspi->annotfile_fn); exit(1);
	}
	if (pspi->datafile_fn[0] && (pspi->datafile = fopen(pspi->datafile_fn, "r")) == 0)
	{
		dcerror(1, "Could not open datafile %s\n", pspi->datafile_fn); exit(1);
	}
	if (pspi->filterfile_fn[0] && (pspi->filterfile = fopen(pspi->filterfile_fn, "r")) == 0)
	{
		dcerror(1, "Could not open filterfile %s\n", pspi->filterfile_fn); exit(1);
	}
	if (pspi->weightfile_fn[0] && pspi->annotfile_fn[0])
		spi->use_func_weights=1;
	else if (pspi->weightfile_fn[0] && !pspi->annotfile_fn[0])
		printf("Warning: weightfile specified but not annotfile so will not assign weights according to function\n");
	else if (!pspi->weightfile_fn[0] && pspi->annotfile_fn[0])
		printf("Warning: annotfile specified but not weightfile so will not assign weights according to function\n");
}

// need to consider that not all subjects may be typed for all loci
// code relies on fact that sub is allocated with calloc, which will have set all alleles to 0
// NB the annot file MUST have two sets of counts - for cases and controls
int read_all_data(par_info*pi,sa_par_info *spi,psa_par_info *pspi,subject **sub,int *nsubptr,char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH],float func_weight[MAX_LOCI])
{
	int nsub,first,s,a,l,func_pos,f;
	char dline[1000],aline[1000],pos[100],rsname[100],ref[100],alls[100],all[2][100],sname[100],ccstr[100],effect[100],*ptr;
	float weight;
	std::map<std::string,int> subIDs;
	std::map<std::string,float> weightMap;
	std::map<std::string,std::string> effectMap;
	if (pspi->weightfile)
	{
		while (fgets(aline, 999, pspi->weightfile))
		{
			sscanf(aline,"%s %f",effect,&weight);
			weightMap[effect]=weight;
		}
	}
	if (pspi->annotfile)
	{
		fgets(aline, 999, pspi->annotfile); // ignore first line (though could use it to see how many cohorts there are)
		for (func_pos=0,ptr=aline;strncmp(ptr,"FUNC",4);ptr=skip_word(ptr))
			if (!*ptr)
			{
				dcerror(1, "Could not find FUNC in annotation file %s:\n%s\n", pspi->annotfile_fn, aline); exit(1);
			}
			else
				++func_pos;
		while (fgets(aline, 999, pspi->annotfile))
		{
			sscanf(aline,"%s",pos);
			ptr=aline;
			for (f=0;f<func_pos;++f)
				ptr=skip_word(ptr);
			sscanf(ptr,"%s",effect);
			effectMap[pos]=effect;
		}
	}
	first=1;
	nsub=0;
	pi->nloci=0;
	while (fgets(dline, 999, pspi->datafile))
	{
		if (sscanf(dline,"%s",pos)!=1)
			continue;
//		if (!strncmp("chr", pos, 3) && strchr(pos, ':')) //looks like new locus, not subject ID, will break if subject does have ID like this
		if (strchr(pos, ':')) //looks like new locus, not subject ID, will break if subject does have ID like this
		{
			++pi->nloci;
			l=pi->nloci-1;
			pi->n_alls[l]=2;
			if (sscanf(dline, "%s %s %s", pos, rsname, alls) != 3)
			{
				dcerror(1, "Could not read this locus line:\n%s\n", dline);
				exit(1);
			}
			if (sscanf(alls,"%[^/]",ref)!=1)
			{
				dcerror(1, "Could not read reference allele in this locus line:\n%s\n", dline);
				exit(1);
			}
			if (pspi->weightfile && pspi->annotfile)
			{
				std::map<std::string,std::string>::const_iterator effIter =effectMap.find(pos);
				if (effIter == effectMap.end())
				{
					dcerror(1,"%s not found in annotation file %s\n",pos,pspi->annotfile_fn);
					exit(1);
				}
				sprintf(aline, "%s_%s", pos, effIter->second.c_str());
				aline[19]='\0';
				strcpy(names[l],aline);
				sprintf(comments[l],"%s_%s_%s_%s", pos, rsname,effIter->second.c_str(),alls);
				std::map<std::string,float>::const_iterator weightIter =weightMap.find(effIter->second);
				if (weightIter==weightMap.end())
					{
						dcerror(1,"weight for effect %s not found in weight file %s\n",effIter->second.c_str(),pspi->weightfile_fn);
						exit(1);
					}
				else 
					func_weight[l]=weightIter->second;
			}
			else
			{
				sprintf(names[l], "%s", pos);
				sprintf(comments[l],"%s_%s", pos, alls);
				func_weight[l]=1;
			}

		}
		else
		{
			if (sscanf(dline,"%s %*d %s %[^/]/%s", sname, ccstr, all[0], all[1]) != 4)
			{
				dcerror(1, "Could not read this genotype line:\n%s\n", dline); exit(1);
			}
			std::map<std::string,int>::const_iterator iter =subIDs.find(sname);
			if (iter == subIDs.end()) // found a subject we haven't seen before
			{
				subIDs[sname]=nsub;
				s=nsub;
				strcpy(sub[s]->id,sname);
				if (!strcmp(ccstr,"CONTROL"))
					sub[s]->cc=0;
				else if (!strcmp(ccstr,"CASE"))
					sub[s]->cc=1;
				else
				{
					dcerror(1, "Unrecognised phenotype in this line:\n%s\n", dline); exit(1);
				}
				++nsub;
			}
			else
				s=iter->second;
			if (all[0][0]=='.')
				continue;
			else for (a=0;a<2;++a)
				if (!strcmp(all[a],ref))
					sub[s]->all[l][a]=1;
				else
					sub[s]->all[l][a]=2;
		}
	}
	pi->n_loci_to_use=pi->nloci;
	for (l=0;l<pi->nloci;++l)
		pi->loci_to_use[l]=l;
	*nsubptr=nsub;
	return 1;
}

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub,l,n_new_sub,real_nsub;
	float *score,p;
	int max_cc[2],s,n_non_mendelian;
	non_mendelian *non_mendelians;
	char *non_mendelian_report;
	par_info pi;
	sa_par_info spi;
	psa_par_info pspi;
	subject **sub,**new_sub,**real_sub;
	printf("%s v%s\n",PROGRAM,PSAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);

	assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
	for (s=0;s<MAX_SUB;++s)
		assert(sub[s]=(subject *)calloc(1,sizeof(subject)));
	assert(score=(float *)calloc(MAX_SUB,sizeof(float)));
	max_cc[0]=max_cc[1]=0;
	make_arg_string(arg_string,argc,argv);
	parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi,&spi,&pspi);
	fprintf(pspi.outfile,"pscoreassoc output\n");
	if (pspi.filterfile)
		initExclusions(pspi.filterfile);
	read_all_data(&pi,&spi,&pspi,sub,&nsub,names,comments,func_weight);

if (spi.use_trios)
{
	if (atoi(comments[0])>22 || toupper(comments[0][0]) == 'X' || toupper(comments[0][0]) == 'Y' ||
		toupper(comments[0][0]) == 'C'&&toupper(comments[0][1]) == 'H'&&toupper(comments[0][2]) == 'R' && (atoi(comments[0] + 3) > 22 || toupper(comments[0][3]) == 'X' || toupper(comments[0][3]) == 'Y'))
	{
		error("Cannot at present use trios for genes on X or Y chromosome", "");
		return 1;
	}
}

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


fprintf(pspi.outfile,
"Locus                   controls     frequency        cases          frequency   frequency allele  weight\n"
"                     AA  :   AB  :  BB                  AA  :   AB  :  BB                     \n");
get_freqs(sub,nsub,&pi,&spi,cc_freq,cc_count,cc_genocount);
applyExclusions(&pi);
set_weights(pspi.outfile,weight,missing_score,rarer,sub,nsub,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
get_scores(score,weight,missing_score,rarer,sub,nsub,&pi);
p=do_score_onetailed_ttest(pspi.outfile,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer);

if (spi.do_recessive_test)
{
if (atoi(comments[0])>22 || toupper(comments[0][0])=='X' || toupper(comments[0][0])=='Y' || 
	toupper(comments[0][0])=='C'&&toupper(comments[0][1])=='H'&&toupper(comments[0][2])=='R'&&(atoi(comments[0]+3)>22 || toupper(comments[0][3])=='X' || toupper(comments[0][3])=='Y'))
// simple trick for now to avoid X and Y genes
	fprintf(pspi.outfile,"\nCannot do recessive test for genes on X or Y chromosome.\n");
else if (spi.use_haplotypes)
	do_recessive_HWE_test_with_haplotypes(pspi.outfile,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
else 
	do_recessive_HWE_test(pspi.outfile,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
}

if (spi.use_trios)
{
	fprintf(pspi.outfile,"\n%s\n",*non_mendelian_report?non_mendelian_report:"No non_mendelian transmissions or de novo mutations found\n");
}

stateExclusions(pspi.outfile);
printf("\nProgram run completed OK\n");
return 0;

}


