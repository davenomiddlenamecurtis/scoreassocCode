#include "scoreassoc.hpp"

void do_recessive_HWE_test_with_haplotypes(FILE *fo, float *score, subject **sub, int nsub, par_info *pi, sa_par_info *spi,float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], float *weight, float *missing, int *old_rarer, char names[MAX_LOCI][20])
{
	float tab[2][2],ex[2][2],col_tot[3],row_tot[2],N,counts[2][3],hom_counts[2],all_count[MAX_LOCI],hapMAF;
	int r,c,s,l,ll,lll,a,first,all_genotype,rec_rarer[MAX_LOCI];
	int n_loci_to_remove,loci_to_remove[MAX_LOCI],has_minor[2],num_minor;
	double p,df,dchi;
	float freq,tot,geno_prob[3],geno_count[MAX_LOCI][3];
	int status,which=1;
	par_info rec_pi; // use this to keep track of which loci we are using
	// first thing we will do is remove loci which fail the weight threshold
	// note that weights and rarer are indexed up to pi->n_loci_to_use
	for (l=0,ll=0;l<pi->n_loci_to_use;++l)
	{
		if (weight[l]>=spi->weight_threshold)
		{
			rec_pi.loci_to_use[ll]=pi->loci_to_use[l];
			rec_rarer[rec_pi.loci_to_use[ll]]=rarer[l]; // here is where we go back to indexing rec_rarer in the same way as other locus attributes
			++ll;
		}
	}
	rec_pi.n_loci_to_use=ll;

	// from now on it impossible to marry up the weights with the loci listed in rec_pi
	// hope that is OK
	for (r=0;r<2;++r)
	{
		for (c=0;c<3;++c)
			counts[r][c]=0;
		row_tot[r]=0;
		hom_counts[r]=0;
	}
	fprintf(fo,"\nRecessive analyses based on loci with weight reaching %f\n",spi->weight_threshold);
	fprintf(fo,"Using the following loci:\n");
	for (l=0;l<rec_pi.n_loci_to_use;++l)
		fprintf(fo,"%3d %s\n",l+1,names[rec_pi.loci_to_use[l]]);
	fprintf(fo,"\nSubjects with minor alleles in different haplotypes:\n");
	num_minor=0; // number of haplotypes which bear a minor allele
	for (s = 0; s < nsub; ++s)
	{
		has_minor[0]=has_minor[1]=0;
        for (l=0;l<rec_pi.n_loci_to_use;++l)
			{
				for (a=0;a<2;++a)
					if (sub[s]->all[rec_pi.loci_to_use[l]][a]==rec_rarer[rec_pi.loci_to_use[l]])
						has_minor[a]=1;
			}
		num_minor=has_minor[0]+has_minor[1];
		++counts[sub[s]->cc][num_minor];
		if (num_minor==2)
		{
			fprintf(fo,"%-10s %d\n",sub[s]->id,sub[s]->cc);
			for (a = 0; a < 2; ++a)
			{
				for (l = 0; l < rec_pi.n_loci_to_use; ++l)
				{
					fprintf(fo,"%d ",1+(sub[s]->all[rec_pi.loci_to_use[l]][a] == rec_rarer[rec_pi.loci_to_use[l]]));
				}
				fprintf(fo,"\n");
			}
		}
	}
	fprintf(fo,"\nCounts of haplotypes bearing a minor allele for cases and controls\n");
	fprintf(fo,"\n(AB have at least one minor allele on one haplotype, BB have minor alleles on both haplotypes)\n");
	fprintf(fo,"            AA      AB      BB   \n");
	for (r=0;r<2;++r)
	{
		fprintf(fo,"%8s ",r==0?"Controls":"Cases");
		for (c=0;c<3;++c)
			fprintf(fo,"%7.2f ",counts[r][c]);
		fprintf(fo,"\n");
	}
	col_tot[0]=col_tot[1]=N=0;
	for (r=0;r<2;++r)
	{
		col_tot[0]+=tab[r][0]=counts[r][0]+counts[r][1];
		col_tot[1]+=tab[r][1]=counts[r][2];
		N+=row_tot[r]=tab[r][0]+tab[r][1];
    }
	fprintf(fo,"%10s","\nObserved:\n");
	for (r=0;r<2;++r)
	{
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",tab[r][c]);
		fprintf(fo,"\n");
	}
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			ex[r][c]=row_tot[r]*col_tot[c]/N;
	fprintf(fo,"%10s","\nExpected:\n");
	for (r=0;r<2;++r)
	{
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",ex[r][c]);
		fprintf(fo,"\n");
	}
	dchi=0;
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			if (ex[r][c]!=0)
				dchi+=(tab[r][c]-ex[r][c])*(tab[r][c]-ex[r][c])/(ex[r][c]<1?1:ex[r][c]); /* if expected less than 1 set it to 1 */
	if (col_tot[0]==0 || col_tot[1]==0)
		dchi=0; 
	p=chistat(dchi,1.0)/2;
	if (ex[1][1]>tab[1][1])
		p=1-p;
	fprintf(fo,"\n\nRecessive chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));

	hapMAF=(counts[0][1]*0.5+counts[0][2]+counts[1][1]*0.5+counts[1][2])/nsub;
	geno_prob[0]=(1-hapMAF)*(1-hapMAF);
	geno_prob[1]=2*hapMAF*(1-hapMAF);
	geno_prob[2]=hapMAF*hapMAF;
	fprintf(fo,"\nExpected probabilities for subjects to have total of 0, 1 or 2 haplotypes bearing a minor allele:\n%8.6f %8.6f %8.6f \n",
	geno_prob[0],geno_prob[1],geno_prob[2]);
	/* here we do the recessive HWE test for all loci combined as compound heterozygotes */
	tab[0][1]=counts[1][2];
	ex[0][1]=geno_prob[2]*row_tot[1];
	tab[0][0]=row_tot[1]-tab[0][1];
	ex[0][0]=row_tot[1]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
			/* using Yates correction */
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\nRecessive HWE chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[1],tab[0][1]-1,geno_prob[2]);
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));
	tab[0][1]=counts[0][2];
	ex[0][1]=geno_prob[2]*row_tot[0];
	tab[0][0]=row_tot[0]-tab[0][1];
	ex[0][0]=row_tot[0]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\n\nRecessive HWE for controls chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[0],tab[0][1]-1,geno_prob[2]);
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE for controls exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));	
}

