#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <map>

#include "scoreassoc.hpp"

#include "dcexpr.hpp"
#include "safilterfuncs.hpp"

#define PROGRAM "pscoreassoc"
#define PSAVERSION "3.2"

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


option opt[]=
{
	{ "psdatafile", PSDATAFILE },
	{ "gcdatafile", GCDATAFILE },
	{ "gendatafile", GENDATAFILE },
	{ "weightfile", WEIGHTFILE },
	{ "filterfile", FILTERFILE },
	{ "annotfile", ANNOTFILE },
	{ "locusfilterfile", LOCUSFILTERFILE },
	{ "locusnamefile", LOCUSNAMEFILE },
	{ "locusweightfile", LOCUSWEIGHTFILE },
	{ "casefreqfile", CASEFREQFILE },
	{ "contfreqfile", CONTFREQFILE },
	{ "samplefile", SAMPLEFILE },
	{ "triofile", TRIOFILE },
	{ "outfile", OUTFILE },
	{ "scorefile", SCOREFILE },
	{ "nostringtomatchthis", NUMDATAFILETYPES },
	{ "numloci", NUMLOCI },
// need this if using gc files
	{ "ldthreshold", LDTHRESHOLD },
	{ "minweight", WEIGHTTHRESHOLD },
	{ "dorecessive", DORECESSIVE },
	{ "usehaps", USEHAPS },
	{ "weightfactor", WEIGHTFACTOR },
	{ "argfile", ARGFILE },
	{"", NUMOPTS}
};
// readable files must be listed before writable files

void usage()
{
	printf("pscoreassoc --psdatafile file || --gcdatafile file || --gendatafile file     [options]\n\nOptions:\n"
"--weightfactor x (weight for very rare variants, default 10)\n"
"--outfile file\n"
"--weightfile file (specify functional weights)\n"
"--locusweightfile file (specify functional weight for each locus)\n"
"--locusnamefile file\n"
"--annotfile file (annotations from plink/seq)\n"
"--filterfile file\n"
"--triofile file\n"
"--locusfilterfile file (exclude specific loci)\n"
"--ldthreshold x (to discard variants in LD for recessive analysis, default 0.9)\n"
"--minweight x (to include in recessive analysis, default 0, i.e. all variants)\n"
"--dorecessive\n"
"--usehaps\n"
"--numloci x (needed with --gcdatafile)\n"
"--argfile file (additional arguments)\n"
);
}

#define MAXDEPTH 5

int getNextArg(char *nextArg, int argc,char *argv[], FILE *fp[MAXDEPTH],int *depth, int *argNum)
{
	*nextArg='\0';
	while (*depth>-1)
	{
		if (fscanf(fp[*depth],"%s ",nextArg)==1)
			return 1;
		else
		{
			fclose(fp[*depth]);
			--*depth;
		}
	}
	if (*argNum < argc)
	{
		strcpy(nextArg,argv[*argNum]);
		++ *argNum;
		return 1;
	}
	else
		return 0;
}

int read_all_args(char *argv[],int argc, par_info *pi, sa_par_info *spi)
{
	char arg[2000];
	int a,arg_depth,arg_num;
	FILE *fp[MAXDEPTH];
	arg_depth=-1;
	arg_num=1;
	pi->nloci=0;
	spi->use_locus_names=spi->use_comments=1;
	spi->wfactor=10;
	spi->do_recessive_test=0;
	spi->LD_threshold=0.9;
	spi->use_haplotypes=spi->use_trios=0;
	spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=0;
	while (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num))
	{
		if (strncmp(arg, "--", 2))
		{
			printf("Need to specify options with --\n%s\n\n", arg);
			usage();
			exit(1);
		}
		for (a=0;a<NUMOPTS-1;++a)
			if (!strncmp(arg+2,opt[a].str,strlen(opt[a].str)))
				break;
		if (a == NUMOPTS - 1)
			{
			printf("Unrecognised option: \n%s\n\n", arg);
			usage();
			exit(1);
			}
		int error=0;
		switch (opt[a].o)
		{
		case ARGFILE:
			if (++arg_depth >= MAXDEPTH)
			{
				dcerror(1, "Attempting to recurse too deeply into arg-files with this one: %s\n", arg);
				return 0;
			}
			else if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || arg[0]=='-')
				error=1;
			else
			{
				fp[arg_depth] = fopen(arg, "r");
				if (fp[arg_depth] == NULL)
				{
					dcerror(1, "Could not open arg file: %s\n", arg);
					return 0;
				}
			}
			break;
		case NUMLOCI:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg, "%d", &pi->nloci) != 1)
				error=1;
			break;
		case WEIGHTFACTOR:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg, "%f", &spi->wfactor) != 1)
				error=1;
			break;
		case LDTHRESHOLD:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg,"%f",&spi->LD_threshold)!=1)
				error=1;
			break;
		case WEIGHTTHRESHOLD:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || sscanf(arg,"%f",&spi->weight_threshold)!=1)
				error=1;
			break;
		case DORECESSIVE:
			spi->do_recessive_test=1;
			break;
		case USEHAPS:
			spi->use_haplotypes=1;
			break;
		case PSDATAFILE:
		case GCDATAFILE:
		case GENDATAFILE:
		case WEIGHTFILE:
		case ANNOTFILE:
		case FILTERFILE:
		case LOCUSFILTERFILE:
		case LOCUSWEIGHTFILE:
		case LOCUSNAMEFILE:
		case SAMPLEFILE:
		case CASEFREQFILE:
		case CONTFREQFILE:
		case OUTFILE:
		case SCOREFILE:
			if (getNextArg(arg, argc, argv, fp,&arg_depth, &arg_num) == 0 || arg[0]=='-' || sscanf(arg, "%s",spi->df[opt[a].o].fn) != 1)
				error=1;
			break;
		}
		if (error)
			{
				dcerror(1, "Error reading %s\n %s", opt[a].str,arg); exit(1);
			}
	}
	return 1;
}

int process_options(par_info *pi, sa_par_info *spi)
{
	int a,l;
	if ((spi->df[PSDATAFILE].fn[0]==0)+(spi->df[GCDATAFILE].fn[0]==0)+(spi->df[GENDATAFILE].fn[0]==0) < 2)
	{
		dcerror(1,"Must specify only one of --psdatafile, --gcdatafile or --gendatafile"); exit(1);
	}
	if (spi->df[GCDATAFILE].fn[0]!='\0' && pi->nloci==0)
	{
		dcerror(1,"Need to specify --nloci when using --gcdatafile"); exit(1);
	}

	for (l = 0; l < pi->nloci; ++l)
		pi->n_alls[l]=2;
	for (a=0;a<OUTFILE;++a)
		if (spi->df[a].fn[0] && (spi->df[a].fp = fopen(spi->df[a].fn, "r")) == 0)
		{
			dcerror(1, "Could not open %s %s\n", opt[a].str,spi->df[a].fn); exit(1);
		}
	for (a=OUTFILE;a<=SCOREFILE;++a)
		if (spi->df[a].fn[0] && (spi->df[a].fp = fopen(spi->df[a].fn, "w")) == 0)
		{
			dcerror(1, "Could not open %s %s\n", opt[a].str,spi->df[a].fn); exit(1);
		}
	if (spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
		spi->use_func_weights=1;
	else if (spi->df[WEIGHTFILE].fp && !spi->df[ANNOTFILE].fp)
		printf("Warning: weightfile specified but not annotfile so will not assign weights according to function\n");
	else if (!spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
		printf("Warning: annotfile specified but not weightfile so will not assign weights according to function\n");
	if (spi->df[LOCUSWEIGHTFILE].fp)
	{
		spi->use_func_weights=1;
		if (spi->df[WEIGHTFILE].fp)
			printf("Warning: weightfile specified but will read weights from locusweightfile instead\n");
	}
	if (spi->df[TRIOFILE].fp)
		spi->use_trios=1;
}

char *skip_word(char *ptr)
{
	while (*ptr && !isspace(*ptr))
		++ptr;
	while (*ptr && isspace(*ptr))
		++ptr;
	return ptr;
}

// need to consider that not all subjects may be typed for all loci
// code relies on fact that sub is allocated with calloc, which will have set all alleles to 0
// NB the annot file MUST have two sets of counts - for cases and controls

int read_ps_datafile(par_info*pi, sa_par_info *spi, subject **sub, int *nsubptr, char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH], float func_weight[MAX_LOCI],	
	std::map<std::string,float> weightMap,	std::map<std::string,std::string> effectMap)
{
	int nsub,first,s,a,l,f;
	char dline[1000],aline[1000],pos[100],rsname[100],ref[100],alls[100],all[2][100],sname[100],ccstr[100],*ptr;
	float weight;
	std::map<std::string,int> subIDs;
	first=1;
	nsub=0;
	pi->nloci=0;
	while (fgets(dline, 999, spi->df[PSDATAFILE].fp))
	{
		if (sscanf(dline,"%s",pos)!=1)
			continue;
		if (!strncmp("chr", pos, 3) && strchr(pos, ':')) //looks like new locus, not subject ID, will break if subject does have ID like this
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
			if (spi->df[WEIGHTFILE].fp && spi->df[ANNOTFILE].fp)
			{
				std::map<std::string,std::string>::const_iterator effIter=effectMap.find(pos);
				if (effIter == effectMap.end())
				{
					dcerror(1,"%s not found in annotation file %s\n",pos,spi->df[ANNOTFILE].fn);
					exit(1);
				}
				sprintf(comments[l],"%s_%s_%s", pos, effIter->second.c_str(),alls);
				std::map<std::string,float>::const_iterator weightIter =weightMap.find(effIter->second);
				if (weightIter==weightMap.end())
					{
						dcerror(1,"weight for effect %s not found in weight file %s\n",effIter->second.c_str(),spi->df[WEIGHTFILE].fn);
						exit(1);
					}
				else 
					func_weight[l]=weightIter->second;
			}
			else
			{
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
	*nsubptr=nsub;
	return 1;
}

int read_all_data(par_info *pi,sa_par_info *spi,subject **sub,int *nsubptr,char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH],float func_weight[MAX_LOCI])
{
	char aline[1000],pos[100],effect[100],*ptr;
	int func_pos,l,f,use;
	float wt;
	std::map<std::string,float> weightMap;
	std::map<std::string,std::string> effectMap;
	if (spi->df[WEIGHTFILE].fp)
	{
		while (fgets(aline, 999, spi->df[WEIGHTFILE].fp))
		{
			sscanf(aline,"%s %f",effect,&wt);
			weightMap[effect]=wt;
		}
	}
	if (spi->df[ANNOTFILE].fp)
	{
		fgets(aline, 999, spi->df[ANNOTFILE].fp); // ignore first line (though could use it to see how many cohorts there are)
		for (func_pos=0,ptr=aline;strncmp(ptr,"FUNC",4);ptr=skip_word(ptr))
			if (!*ptr)
			{
				dcerror(1, "Could not find FUNC in annotation file %s:\n%s\n", spi->df[ANNOTFILE].fn, aline); exit(1);
			}
			else
				++func_pos;
		while (fgets(aline, 999, spi->df[ANNOTFILE].fp))
		{
			sscanf(aline,"%s",pos);
			ptr=aline;
			for (f=0;f<func_pos;++f)
				ptr=skip_word(ptr);
			sscanf(ptr,"%s",effect);
			effectMap[pos]=effect;
		}
	}
	if (spi->df[PSDATAFILE].fp)
		read_ps_datafile(pi,spi,sub,nsubptr,names,comments, func_weight,weightMap,effectMap);
	else if (spi->df[GCDATAFILE].fp)
		read_all_subjects(spi->df[GCDATAFILE].fp,sub,nsubptr,pi);
	if (spi->df[LOCUSFILTERFILE].fp)
	{
		pi->n_loci_to_use=0;
		for (l = 0; l < pi->nloci; ++l)
		{
			if (fscanf(spi->df[LOCUSFILTERFILE].fp, " %d", &use) != 1)
			{
				dcerror(1, "Not enough values in locusfilterfile %s\n", spi->df[LOCUSFILTERFILE].fn); exit(1);
			}
			else if (use==1)
				pi->loci_to_use[pi->n_loci_to_use++]=l;
		}
	}
	else
	{
		pi->n_loci_to_use = pi->nloci;
		for (l = 0; l < pi->nloci; ++l)
			pi->loci_to_use[l] = l;
	}
	if (spi->df[LOCUSWEIGHTFILE].fp)
	{
		for (l = 0; l < pi->nloci; ++l)
			if (fscanf(spi->df[LOCUSWEIGHTFILE].fp,"%f ",&func_weight[l])!=1)
			{
				dcerror(1, "Not enough values in locusweightfile %s\n", spi->df[LOCUSWEIGHTFILE].fn); exit(1);
			}
	}
	else if (spi->df[WEIGHTFILE].fp ==0)
		for (l = 0; l < pi->nloci; ++l)
			func_weight[l]=1;
	if (spi->df[LOCUSNAMEFILE].fp)
	{
		for (l = 0; l < pi->nloci; ++l)
			if (fscanf(spi->df[LOCUSNAMEFILE].fp,"%s ",&comments[l])!=1)
			{
				dcerror(1, "Not enough values in locusnamefile %s\n", spi->df[LOCUSNAMEFILE].fn); exit(1);
			}
	}
	for (l = 0; l < pi->nloci; ++l)
	{
		strncpy(names[l],comments[l],NAME_LENGTH-1);
		names[l][NAME_LENGTH-1]='\0';
	}

	return 1;
}

int main(int argc, char *argv[])
{
	char arg_string[2000];
	int nsub,n_new_sub,real_nsub;
	float *score,p;
	int max_cc[2],s,n_non_mendelian;
	non_mendelian *non_mendelians;
	char *non_mendelian_report;
	par_info pi;
	sa_par_info spi;
	subject **sub,**new_sub,**real_sub;
	pi.use_cc=1;
	printf("%s v%s\n",PROGRAM,PSAVERSION);
	printf("MAX_LOCI=%d\nMAX_SUB=%d\n",MAX_LOCI,MAX_SUB);

	assert(sub=(subject **)calloc(MAX_SUB,sizeof(subject*)));
	for (s=0;s<MAX_SUB;++s)
		assert(sub[s]=(subject *)calloc(1,sizeof(subject)));
	assert(score=(float *)calloc(MAX_SUB,sizeof(float)));
	max_cc[0]=max_cc[1]=0;
	read_all_args(argv,argc, &pi, &spi);
	// make_arg_string(arg_string,argc,argv);
	// parse_arg_string(arg_string,&pi,&spi,&pspi);
	process_options(&pi,&spi);
	if (spi.df[FILTERFILE].fp)
		initExclusions(spi.df[FILTERFILE].fp);
	read_all_data(&pi,&spi,sub,&nsub,names,comments,func_weight);
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
	fprintf(spi.df[OUTFILE].fp,"pscoreassoc output\n"
"Locus                   controls     frequency        cases          frequency   frequency allele  weight\n"
"                     AA  :   AB  :  BB                  AA  :   AB  :  BB                     \n");
get_freqs(sub,nsub,&pi,&spi,cc_freq,cc_count,cc_genocount);
applyExclusions(&pi);
set_weights(spi.df[OUTFILE].fp,weight,missing_score,rarer,sub,nsub,&pi,&spi,func_weight,cc_freq,cc_count,max_cc,names,comments);
get_scores(score,weight,missing_score,rarer,sub,nsub,&pi);
p=do_score_onetailed_ttest(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer);
if (spi.df[SCOREFILE].fp)
	write_scores(spi.df[SCOREFILE].fp,sub,nsub,score);
if (spi.do_recessive_test)
{
if (atoi(comments[0]+3)>22 || toupper(comments[0][4])=='X' || toupper(comments[0][4])=='Y')
// simple trick for now to avoid X and Y genes
	fprintf(spi.df[OUTFILE].fp,"\nCannot do recessive test for genes on X or Y chromosome.\n");
else if (spi.use_haplotypes)
	do_recessive_HWE_test_with_haplotypes(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
else 
	do_recessive_HWE_test(spi.df[OUTFILE].fp,score,sub,nsub,&pi,&spi,cc_freq,cc_count,max_cc,weight,missing_score,rarer,names);
}

stateExclusions(spi.df[OUTFILE].fp);
printf("\nProgram run completed OK\n");
return 0;

}


