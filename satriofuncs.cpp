#include "scoreassoc.hpp"


int sort_trios(subject **sub, int nsub, par_info *pi, subject **new_sub, non_mendelian *non_mendelians,int *n_non_mendelian,char *non_mendelian_report)
{
	FILE *ft,*flog;
	char line[201],f_id[MAX_ID_LENGTH],c_id[MAX_ID_LENGTH],m_id[MAX_ID_LENGTH];
	subject *child_ptr,*parent_ptr[2],*case_ptr,*cont_ptr;
	int s,n_new_sub,pl,l,d,p,a,n_in_parents,n_in_this_parent,has_allele[2],test_allele,assigned_parent[2],is_non_mendelian[2],child_num,is_de_novo;
	*non_mendelian_report='\0';
	*n_non_mendelian=0;
	ft=fopen(trios_fn,"r");
	if (ft == NULL)
	{
		error("Could not open file giving pedigree information about trios: ",trios_fn);
		return 0;
	}
	flog=fopen("trios.log","w");
	n_new_sub=0;
	while (fgets(line, 200, ft))
	{
		if (sscanf(line,"%s %s %s",c_id,f_id,m_id)!=3)
			break;
		if (!strcmp(f_id, "0"))
		{
			if (strcmp(m_id, "0"))
			{
				error("Paternal ID is 0 but maternal ID is not, in line: ", line);
				return 0;
			}
			else
				continue;
		}
		child_ptr=parent_ptr[0]=parent_ptr[1]=0;
		for (s = 0; s < nsub; ++s)
		{
			if (!child_ptr && !strcmp(sub[s]->id, c_id))
			{
				child_ptr = sub[s];
				child_num = s;
			}

			if (!parent_ptr[0] && !strcmp(sub[s]->id,f_id))
				parent_ptr[0]=sub[s];
			if (!parent_ptr[1] && !strcmp(sub[s]->id,m_id))
				parent_ptr[1]=sub[s];
		}
		if (!child_ptr)
		{
			error("When searching information from pedigree/trios file, could not find subject with this ID in data file: ",c_id);
			return 0;
		}
		if (!parent_ptr[0])
		{
			error("When searching information from pedigree/trios file, could not find subject with this ID in data file: ",f_id);
			return 0;
		}
		if (!parent_ptr[1])
		{
			error("When searching information from pedigree/trios file, could not find subject with this ID in data file: ",m_id);
			return 0;
		}
		case_ptr=new_sub[n_new_sub];
		cont_ptr=new_sub[n_new_sub+1];
		sprintf(case_ptr->id,"%s_CASE",child_ptr->id);
		sprintf(cont_ptr->id,"%s_CONT",child_ptr->id);
		case_ptr->cc=1;
		cont_ptr->cc=0;
		for (pl = 0; pl < pi->n_loci_to_use; ++pl)
		{
			l=pi->loci_to_use[pl];
			for (a=0;a<2;++a)
				case_ptr->all[l][a]=cont_ptr->all[l][a]=0;
			if (parent_ptr[0]->all[l][0] == 0 || parent_ptr[1]->all[l][0] == 0 || child_ptr->all[l][0] == 0)
			{
				continue;
				// assume all in trio unknown if one unknown - may anyway indicate a problem genotyping this variant
			}
			for (a = 0; a < 2; ++a)
				is_non_mendelian[a]=assigned_parent[a]=0;
			is_de_novo=0;
			// first check for typical de novo, i.e. where both parents 11
			for (p = 0; p < 2; ++p)
			{
				for (a=0;a<2;++a)
					if (parent_ptr[p]->all[l][a]!=1)
						goto not_de_novo;
			}
			if (child_ptr->all[l][0] == 1 && child_ptr->all[l][1] == 1)
			{
				case_ptr->all[l][0]=case_ptr->all[l][1]=cont_ptr->all[l][0]=cont_ptr->all[l][1]=1;
				goto sorted_trio;
			}
			else
			{
				cont_ptr->all[l][0]=cont_ptr->all[l][1]=1;
				for (a = 0; a < 2; ++a)
					case_ptr->all[l][a] = child_ptr->all[l][a];
				is_de_novo=1;
			}

			not_de_novo:
			// deal with homozygous child next
			if ((test_allele = child_ptr->all[l][0]) == child_ptr->all[l][1])
			{
				for (p = 0; p < 2; ++p)
				{
					for (a = 0; a < 2; ++a)
						if (parent_ptr[p]->all[l][a] == test_allele)
						{
							case_ptr->all[l][p]=test_allele;
							cont_ptr->all[l][p]=parent_ptr[p]->all[l][(a+1)%2];
							assigned_parent[p]=1;
							break;
						}
					if (a == 2) // non-mendelian
					{
						is_non_mendelian[p]=1;
						case_ptr->all[l][p]=test_allele;
						cont_ptr->all[l][p]=parent_ptr[p]->all[l][0];
						assigned_parent[p]=1;
					}
				}
			goto sorted_trio;
			}
			// does each allele definitely come from one parent?
			for (d = 0; d < 2; ++d)
			{
				test_allele=child_ptr->all[l][d];
				n_in_parents=0;
				for (a = 0; a < 2; ++a)
					has_allele[a]=0;
				for (p = 0; p < 2; ++p)
				{
					if (assigned_parent[p])
						continue;
					n_in_this_parent = 0;
					for (a = 0; a < 2; ++a)
					{
						if (test_allele == parent_ptr[p]->all[l][a])
						{
							++n_in_parents;
							++n_in_this_parent;
							has_allele[p] = 1;
						}
					}
					if (n_in_this_parent == 2)
					{
						// homozygous parent - assign the allele to this parent and skip other tests
						case_ptr->all[l][p] = cont_ptr->all[l][p] = test_allele;
						assigned_parent[p] = 1;
						goto sorted_one_parent;
					}
				}
				if (n_in_parents == 1) // unambiguous
				{
					for (p=0;p<2;++p)
						if (has_allele[p]) // find out which parent had it
						{
						case_ptr->all[l][p]=test_allele;
						for (a=0;a<2;++a)
							if (parent_ptr[p]->all[l][a] == test_allele)
							{
							cont_ptr->all[l][p]=parent_ptr[p]->all[l][(a+1)%2];
							break;
							}
						assigned_parent[p]=1;
						goto sorted_one_parent;
						}
				}
				else if (n_in_parents == 0) // de novo
				{
					is_non_mendelian[d]=1;
					// assign later when have had a chance to check out other allele to find out which parent it is from
				}
				sorted_one_parent:
					;
				
			}
			// now I have looked at both alleles in the child and tried to assign them to parent
			if (assigned_parent[0] == 0 || assigned_parent[1] == 0) // more to do
			{
				if (is_non_mendelian[0] == 1 && is_non_mendelian[1] == 1) // three homozygotes, e.g. 11, 22, 22
				{
					cont_ptr->all[l][0]=cont_ptr->all[l][1]=parent_ptr[0]->all[l][0];
					case_ptr->all[l][0]=case_ptr->all[l][1]=(parent_ptr[0]->all[l][0]+1)%2;
				}
				else if (is_non_mendelian[0] == 1 || is_non_mendelian[1] == 1) // one de novo, one parent will be sorted
				{
					for (d = 0; d < 2; ++d)
					{
						if (is_non_mendelian[d])
						{
							for (p=0;p<2;++p)
								if (assigned_parent[p] == 0)
								{
								case_ptr->all[l][p]=parent_ptr[p]->all[l][0]%2+1;
								cont_ptr->all[l][p]=parent_ptr[p]->all[l][0];
								}
						}
					}
				}
				else if (assigned_parent[0] ==0 && assigned_parent[1] == 0) // must be 12, 12, 12
				{
					test_allele=case_ptr->all[l][0]=cont_ptr->all[l][1]=parent_ptr[0]->all[l][0];
					case_ptr->all[l][1]=cont_ptr->all[l][0]=test_allele%2+1;
				}
				else if (assigned_parent[0] == 0) // happens if first parent homozygous, child and second heterozygous
				{
					if ((test_allele = parent_ptr[0]->all[l][0]) != parent_ptr[0]->all[l][1])
					{
						error("Bad logic sorting out trio genotypes in trio containing: ", child_ptr->id);
						return 0;
					}
					else
					{
						case_ptr->all[l][0]=cont_ptr->all[l][0]=test_allele;
					}
				}
				// I think that is everything taken care of
				else
				{
					error("Bad logic sorting out trio genotypes in trio containing: ", child_ptr->id);
					return 0;
				}
			}
		sorted_trio:
			if (child_ptr->all[l][0]+child_ptr->all[l][1] != case_ptr->all[l][0]+case_ptr->all[l][1])
			{
				error("Made an error in sorting out case genotype in trio containing: ", child_ptr->id);
				return 0;
			}
			if (parent_ptr[0]->all[l][0]+parent_ptr[0]->all[l][1]+parent_ptr[1]->all[l][0]+parent_ptr[1]->all[l][1] != 
				case_ptr->all[l][0]+case_ptr->all[l][1]+cont_ptr->all[l][0]+cont_ptr->all[l][1]
				&& is_non_mendelian[0]==0 && is_non_mendelian[1]==0 && is_de_novo==0)
			{
				error("Made an error in sorting out control genotype in trio containing: ", child_ptr->id);
				return 0;
			}
			if (is_de_novo)
			{
				sprintf(strchr(non_mendelian_report, '\0'),"De novo mutation detected in transmissions from %s and %s to %s, genotypes: %s:%d%d %s:%d%d %s:%d%d  %s\n",
					parent_ptr[0]->id, parent_ptr[1]->id, child_ptr->id,
					child_ptr->id,child_ptr->all[l][0],child_ptr->all[l][1],
					parent_ptr[0]->id,parent_ptr[0]->all[l][0],parent_ptr[0]->all[l][1],
					parent_ptr[1]->id,parent_ptr[1]->all[l][0],parent_ptr[1]->all[l][1],
					names[l]);
				non_mendelians[*n_non_mendelian].sub=child_num;
				non_mendelians[*n_non_mendelian].loc=l;
				non_mendelians[*n_non_mendelian].nd=DE_NOVO;
				++*n_non_mendelian;
			}
			else for (p = 0; p < 2; ++p)
			{
				if (is_non_mendelian[p])
				{
					sprintf(strchr(non_mendelian_report, '\0'),"Non-mendelian transmission from %s to %s, genotypes: %s:%d%d %s:%d%d %s:%d%d  %s\n",
						parent_ptr[p]->id, child_ptr->id,
						child_ptr->id,child_ptr->all[l][0],child_ptr->all[l][1],
						parent_ptr[0]->id,parent_ptr[0]->all[l][0],parent_ptr[0]->all[l][1],
						parent_ptr[1]->id,parent_ptr[1]->all[l][0],parent_ptr[1]->all[l][1],
						names[l]);
					non_mendelians[*n_non_mendelian].sub=child_num;
					non_mendelians[*n_non_mendelian].loc=l;
					non_mendelians[*n_non_mendelian].nd=NON_MENDELIAN;
					++*n_non_mendelian;
				}
			}
		}
		n_new_sub+=2;
		fprintf(flog,"%d%d %d%d %d%d => %d%d %d%d \n",
			child_ptr->all[l][0],child_ptr->all[l][1],
			parent_ptr[0]->all[l][0],parent_ptr[0]->all[l][1],
			parent_ptr[1]->all[l][0],parent_ptr[1]->all[l][1],
			case_ptr->all[l][0],case_ptr->all[l][1],
			cont_ptr->all[l][0],cont_ptr->all[l][1]);
	}
	fclose(ft);
	fclose(flog);
	return n_new_sub;
}
