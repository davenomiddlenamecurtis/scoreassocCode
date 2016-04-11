#ifndef SCOREASSOCHPP
#define SCOREASSOCHPP

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

extern "C" {
#include "sagcutils.h"
};

#define MAX_COMMENT_LENGTH 100
#define NAME_LENGTH 20

enum OPT {
	PSDATAFILE=0,GCDATAFILE,GENDATAFILE,WEIGHTFILE,ANNOTFILE,FILTERFILE,LOCUSFILTERFILE,LOCUSNAMEFILE,LOCUSWEIGHTFILE,TRIOFILE,SAMPLEFILE,CASEFREQFILE,CONTFREQFILE,OUTFILE,SCOREFILE,NUMDATAFILETYPES,NUMLOCI,LDTHRESHOLD,WEIGHTTHRESHOLD,DORECESSIVE,USEHAPS,WEIGHTFACTOR,ARGFILE,NUMOPTS
};

struct option_t {
	char *str;
	OPT o;
};

typedef struct option_t option;
extern option opt[];

class sa_data_file_type {
public:
	FILE *fp;
	char fn[200];
	sa_data_file_type() { fp = 0; fn[0] = '\0'; }
	~sa_data_file_type() { if (fp != 0 && fp != stdout) fclose(fp); }
};

struct sa_par_info_t {
sa_data_file_type df[NUMDATAFILETYPES];
float wfactor,LD_threshold,weight_threshold;
int use_func_weights,use_cc_freqs[2],use_locus_names,use_comments,do_recessive_test,use_haplotypes,use_trios,use_probs;
};

typedef struct sa_par_info_t sa_par_info;

enum { DE_NOVO=0, NON_MENDELIAN };
struct non_mendelian_t {
	int loc,sub,nd;
};
typedef struct non_mendelian_t non_mendelian;

float do_score_onetailed_ttest(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *rarer);
void get_scores(float *score,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi);
void set_weights(FILE *fo,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
int read_score_assoc_par(FILE *fp,par_info *pi,float *wfactor,int *use_func_weights,int use_cc_freqs[2],float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],int *use_locus_names,char names[MAX_LOCI][20],int *use_comments,char comments[MAX_LOCI][MAX_COMMENT_LENGTH], int *do_recessive_test,float *weight_threshold,float *LD_threshold,int *use_haplotypes);
int read_sa_par(FILE *fp,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_gencount[2][3][MAX_LOCI]);
void write_scores(FILE *fs,subject **sub,int nsub,float *score);
void do_recessive_HWE_test(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *old_rarer,char names[MAX_LOCI][20]);
void do_recessive_HWE_test_with_haplotypes(FILE *fo, float *score, subject **sub, int nsub, par_info *pi,sa_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], float *weight, float *missing, int *old_rarer, char names[MAX_LOCI][20]);
extern int sort_trios(subject **sub,int nsub,par_info *pi,subject **new_sub,non_mendelian *nm,int *n_non_mendelian,char *non_mendelian_report);
int read_all_gen_subjects(FILE *fi,subject **s,int *nsub,par_info *pi);

extern double cumulBinom(int N,int k,double p);

extern float weight[MAX_LOCI],missing_score[MAX_LOCI],func_weight[MAX_LOCI],cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
extern int rarer[MAX_LOCI];
extern char names[MAX_LOCI][NAME_LENGTH],comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];

#define USEFILTERS 1
#endif