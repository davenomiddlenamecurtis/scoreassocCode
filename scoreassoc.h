#ifndef SCOREASSOCH
#define SCOREASSOCH

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "gcutils.h"

#define MAX_COMMENT_LENGTH 100

struct sa_par_info_t {
float wfactor,LD_threshold,weight_threshold;
int use_func_weights,use_cc_freqs[2],use_locus_names,use_comments,do_recessive_test,use_haplotypes,use_trios;
};

typedef struct sa_par_info_t sa_par_info;

struct non_mendelian_t {
	int loc,sub;
	enum { DE_NOVO, NON_MENDELIAN } nd;
};
typedef struct non_mendelian_t non_mendelian;

double do_score_onetailed_ttest(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *rarer);
void get_scores(float *score,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi);
void set_weights(FILE *fo,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
int read_score_assoc_par(FILE *fp,par_info *pi,float *wfactor,int *use_func_weights,int use_cc_freqs[2],float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],int *use_locus_names,char names[MAX_LOCI][20],int *use_comments,char comments[MAX_LOCI][MAX_COMMENT_LENGTH], int *do_recessive_test,float *weight_threshold,float *LD_threshold,int *use_haplotypes);
int read_sa_par(FILE *fp,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH]);
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_gencount[2][3][MAX_LOCI]);
void write_scores(char *fn,subject **sub,int nsub, float *score);
void do_recessive_HWE_test(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *old_rarer,char names[MAX_LOCI][20]);
void do_recessive_HWE_test_with_haplotypes(FILE *fo, float *score, subject **sub, int nsub, par_info *pi, sa_par_info *spi, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], float *weight, float *missing, int *old_rarer, char names[MAX_LOCI][20]);
double do_score_onetailed_wilcoxon(FILE *fo, float *score, subject **sub, int nsub);
extern int sort_trios(subject **sub,int nsub,par_info *pi,subject **new_sub,non_mendelian *nm,int *n_non_mendelian,char *non_mendelian_report);

extern double cumulBinom(int N,int k,double p);
extern double one_tailed_binomial_p(int N, int k, double p);

extern float weight[MAX_LOCI],missing_score[MAX_LOCI],func_weight[MAX_LOCI],cc_freq[2][MAX_LOCI],cc_count[2][MAX_LOCI],cc_genocount[2][3][MAX_LOCI];
extern int rarer[MAX_LOCI];
extern char names[MAX_LOCI][20],comments[MAX_LOCI][MAX_COMMENT_LENGTH],trios_fn[500];

#define USEFILTERS 1
#endif

