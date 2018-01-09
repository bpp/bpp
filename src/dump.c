/*
    Copyright (C) 2016 Tomas Flouri and Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "bpp.h"

#define PROP_COUNT 5

static BYTE dummy[256] = {0};

static void dump_chk_header(FILE * fp, stree_t * stree)
{
  long i;

  BYTE header[76];
  BYTE size_type;

  /* write magic */
  header[0] = BPP_MAGIC[0];
  header[1] = BPP_MAGIC[1];
  header[2] = BPP_MAGIC[2];
  header[3] = BPP_MAGIC[3];

  /* write version major */
  header[4] = (BYTE)((VERSION_MAJOR >>  0) & 0xFF);
  header[5] = (BYTE)((VERSION_MAJOR >>  8) & 0xFF);
  header[6] = (BYTE)((VERSION_MAJOR >> 16) & 0xFF);
  header[7] = (BYTE)((VERSION_MAJOR >> 24) & 0xFF);
  
  /* write version minor */
  header[8]  = (BYTE)((VERSION_MINOR >>  0) & 0xFF);
  header[9]  = (BYTE)((VERSION_MINOR >>  8) & 0xFF);
  header[10] = (BYTE)((VERSION_MINOR >> 16) & 0xFF);
  header[11] = (BYTE)((VERSION_MINOR >> 24) & 0xFF);

  /* write version patch */
  header[12] = (BYTE)((VERSION_PATCH >>  0) & 0xFF);
  header[13] = (BYTE)((VERSION_PATCH >>  8) & 0xFF);
  header[14] = (BYTE)((VERSION_PATCH >> 16) & 0xFF);
  header[15] = (BYTE)((VERSION_PATCH >> 24) & 0xFF);

  /* write checkpoint version */
  header[16] = (BYTE)((VERSION_CHKP >>  0) & 0xFF);
  header[17] = (BYTE)((VERSION_CHKP >>  8) & 0xFF);
  header[18] = (BYTE)((VERSION_CHKP >> 16) & 0xFF);
  header[19] = (BYTE)((VERSION_CHKP >> 24) & 0xFF);

  fwrite((void *)header, sizeof(BYTE), 20, fp);

  size_t size_int = sizeof(int);
  size_t size_long = sizeof(long);
  size_t size_double = sizeof(double);

  /* write size of int */
  size_type = (BYTE)(size_int & 0xFF);
  fwrite((void *)&size_type,sizeof(BYTE),1,fp);

  /* write size of long */
  size_type = (BYTE)(size_long & 0xFF);
  fwrite((void *)&size_type,sizeof(BYTE),1,fp);

  /* write size of double */
  size_type = (BYTE)(size_double & 0xFF);
  fwrite((void *)&size_type,sizeof(BYTE),1,fp);

  /* write RNG value */
  unsigned int rng = get_legacy_rndu_status();
  fwrite((void *)&rng, sizeof(unsigned int), 1, fp);

  /* number of sections */
  unsigned int sections = 3;
  fwrite((void *)&sections, sizeof(unsigned int), 1, fp);

  /* compute length of section 1 */
  unsigned long size_section = 0;
  size_section += sizeof(unsigned long);              /* seed */
  size_section += strlen(opt_msafile)+1;              /* msa filename */
  size_section += strlen(opt_mapfile)+1;              /* imap filename */
  size_section += strlen(opt_outfile)+1;              /* output filename */
  size_section += strlen(opt_mcmcfile)+1;             /* mcmc filename */
  
  size_section += 2*sizeof(long) + 2*sizeof(double);  /* speciesdelimitation */

  size_section += sizeof(long) + 3*sizeof(double);    /* speciestree */

  size_section += sizeof(long);                       /* speciesmodelprior */

  size_section += sizeof(long);                       /* species&tree */
  for (i = 0; i < stree->tip_count; ++i)
    size_section += strlen(stree->nodes[i]->label)+1;

  size_section += sizeof(long);                       /* usedata */
  size_section += sizeof(long);                       /* cleandata */
  size_section += sizeof(long);                       /* nloci */

  size_section += 2*sizeof(double) + sizeof(long);    /* thetaprior */
  size_section += 2*sizeof(double);                   /* tauprior */

  size_section += sizeof(long) + 5*sizeof(double);    /* finetune */
  size_section += stree->tip_count*sizeof(long);      /* diploid */
  size_section += sizeof(long);                       /* burnin */
  size_section += sizeof(long);                       /* sampfreq */
  size_section += sizeof(long);                       /* nsample */
  size_section += sizeof(long);                       /* current step */
  size_section += sizeof(long);                       /* finetune round */
  size_section += sizeof(unsigned long);              /* MCMC file offset */

  size_section += PROP_COUNT*sizeof(double);          /* pjump */

  /* write section 1 size */
  fwrite((void *)&size_section,sizeof(unsigned long),1,fp);
}

static void dump_chk_section_1(FILE * fp, stree_t * stree, double * pjump, long curstep, long ft_round, long mcmc_offset)
{
  size_t i;
//  size_t long_size = sizeof(long);
//  size_t double_size = sizeof(double);


  /* write seed */
  fwrite((void *)&opt_seed,sizeof(long),1,fp);

  /* write seqfile */
  fwrite((void *)opt_msafile, sizeof(BYTE), strlen(opt_msafile)+1, fp);

  /* write imap file */
  fwrite((void *)opt_mapfile, sizeof(BYTE), strlen(opt_mapfile)+1, fp);

  /* write outfile */
  fwrite((void *)opt_outfile, sizeof(BYTE), strlen(opt_outfile)+1, fp);

  /* write mcmcfile */
  fwrite((void *)opt_mcmcfile, sizeof(BYTE), strlen(opt_mcmcfile)+1, fp);

  /* write speciesdelimitation */
  fwrite((void *)&opt_delimit,sizeof(long),1,fp);
  fwrite((void *)&opt_rjmcmc_method,sizeof(long),1,fp);

  /* always write two more elements, even if speciesdelimitation = 0 */
  if (opt_rjmcmc_method == 0)
  {
    fwrite((void *)&opt_rjmcmc_epsilon,sizeof(double),1,fp);

    /* dummy bytes to fill the file */
    fwrite((void *)dummy,sizeof(double),1,fp);
  }
  else
  {
    fwrite((void *)&opt_rjmcmc_alpha,sizeof(double),1,fp);
    fwrite((void *)&opt_rjmcmc_mean,sizeof(double),1,fp);
  }

  /* write speciestree */
  fwrite((void *)&opt_stree,sizeof(long),1,fp);
  fwrite((void *)dummy,3*sizeof(double),1,fp);
  
  /* write speciesmodelprior */
  fwrite((void *)&opt_delimit_prior,sizeof(long),1,fp);

  /* write species&tree */
  fwrite((void *)&(stree->tip_count),sizeof(unsigned int),1,fp);
  for (i = 0; i < stree->tip_count; ++i)
    fwrite((void *)(stree->nodes[i]->label),
           sizeof(char),
           strlen(stree->nodes[i]->label)+1,
           fp);
  
  /* write usedata, cleandata and nloci */
  fwrite((void *)&opt_usedata,sizeof(long),1,fp);
  fwrite((void *)&opt_cleandata,sizeof(long),1,fp);
  fwrite((void *)&opt_nloci,sizeof(long),1,fp);

  /* write theta prior */
  fwrite((void *)&opt_theta_alpha,sizeof(double),1,fp);
  fwrite((void *)&opt_theta_beta,sizeof(double),1,fp);
  fwrite((void *)&opt_est_theta,sizeof(long),1,fp);

  /* write tau prior */
  fwrite((void *)&opt_tau_alpha,sizeof(double),1,fp);
  fwrite((void *)&opt_tau_beta,sizeof(double),1,fp);

  /* write finetune */
  fwrite((void *)&opt_finetune_reset,sizeof(long),1,fp);
  fwrite((void *)&opt_finetune_gtage,sizeof(double),1,fp);
  fwrite((void *)&opt_finetune_gtspr,sizeof(double),1,fp);
  fwrite((void *)&opt_finetune_theta,sizeof(double),1,fp);
  fwrite((void *)&opt_finetune_tau,sizeof(double),1,fp);
  fwrite((void *)&opt_finetune_mix,sizeof(double),1,fp);

  /* write diploid */
  if (opt_diploid)
    fwrite((void *)opt_diploid,sizeof(long),stree->tip_count,fp);
  else
  {
    for (i = 0; i < stree->tip_count; ++i)
      fwrite((void *)dummy,sizeof(long),1,fp);
  }

  /* write mcmc run info */
  fwrite((void *)&opt_burnin,sizeof(long),1,fp);
  fwrite((void *)&opt_samplefreq,sizeof(long),1,fp);
  fwrite((void *)&opt_samples,sizeof(long),1,fp);
  fwrite((void *)&curstep,sizeof(long),1,fp);
  fwrite((void *)&ft_round,sizeof(long),1,fp);

  /* write pjump */
  fwrite((void *)pjump,sizeof(double),PROP_COUNT,fp);

  /* write MCMC file offset */
  fwrite((void *)&mcmc_offset,sizeof(long),1,fp);
}


static void dump_chk_section_2(FILE * fp, stree_t * stree)
{
  long i,j;

  /* write left child node indices */
  for (i = 0; i < stree->inner_count; ++i)
    fwrite((void *)&(stree->nodes[stree->tip_count+i]->left->node_index),
           sizeof(unsigned int),
           1,
           fp);

  /* write right child node indices */
  for (i = 0; i < stree->inner_count; ++i)
    fwrite((void *)&(stree->nodes[stree->tip_count+i]->right->node_index),
           sizeof(unsigned int),
           1,
           fp);

//  /* write branch lengths - TODO: Candidate for removal */
//  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
//    fwrite((void *)&(stree->nodes[i]->length),sizeof(double),1,fp);

  /* write theta */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    fwrite((void *)&(stree->nodes[i]->theta),sizeof(double),1,fp);

  /* write tau */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    fwrite((void *)&(stree->nodes[i]->tau),sizeof(double),1,fp);

  /* write support */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    fwrite((void *)&(stree->nodes[i]->support),sizeof(double),1,fp);

  /* write number of coalescent events */
  assert(opt_nloci == stree->locus_count);
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    fwrite((void *)(stree->nodes[i]->event_count),sizeof(int),opt_nloci,fp);

//  /* write MSC density contribution - TODO: Candidate for removal */
//  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
//    fwrite((void *)&(stree->nodes[i]->logpr_contrib),sizeof(double),1,fp);

  fwrite((void *)&(stree->root_age),sizeof(double),1,fp);

  /* TODO : Perhaps write only seqin_count for tips? */
  /* write number of incoming sequences for each node */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    fwrite((void *)(stree->nodes[i]->seqin_count),sizeof(int),opt_nloci,fp);

  /* write event indices for each node */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    for (j = 0; j < opt_nloci; ++j)
    {
      dlist_item_t * di = stree->nodes[i]->event[j]->head;
      while (di)
      {
        gnode_t * gnode = (gnode_t *)(di->data);

        fwrite((void *)&(gnode->node_index),sizeof(unsigned int),1,fp);

        di = di->next;
      }
    }

  }
}

static void dump_gene_tree(FILE * fp, gtree_t * gtree)
{
  long i;

  /* write gene tree tip labels */
  for (i = 0; i < gtree->tip_count; ++i)
    fwrite((void *)(gtree->nodes[i]->label),
           sizeof(char),
           strlen(gtree->nodes[i]->label)+1,
           fp);

  /* write left child node indices */
  for (i = 0; i < gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[gtree->tip_count+i]->left->node_index),
           sizeof(unsigned int),
           1,
           fp);

  /* write right child node indices */
  for (i = 0; i < gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[gtree->tip_count+i]->right->node_index),
           sizeof(unsigned int),
           1,
           fp);

  /* write branch lengths - TODO: Candidate for removal */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->length),sizeof(double),1,fp);

  /* write ages */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->time),sizeof(double),1,fp);

  /* write population index (corresponding species tree node index) */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->pop->node_index),sizeof(unsigned int),1,fp);

  /* write CLV indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->clv_index),sizeof(unsigned int),1,fp);

  /* write scaler indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->scaler_index),sizeof(int),1,fp);

  /* write pmatrix indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->pmatrix_index),sizeof(unsigned int),1,fp);

  /* write mark - TODO: Candidate for removal */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    fwrite((void *)&(gtree->nodes[i]->mark),sizeof(int),1,fp);

}

static void dump_locus(FILE * fp, gtree_t * gtree, locus_t * locus)
{

  long i;

  /* write number of sites */
  fwrite((void *)&(locus->sites),sizeof(unsigned int),1,fp);

  /* write number of states */
  fwrite((void *)&(locus->states),sizeof(unsigned int),1,fp);

  /* write number of rate categories */
  fwrite((void *)&(locus->rate_cats),sizeof(unsigned int),1,fp);

  /* write number of rate matrices */
  fwrite((void *)&(locus->rate_matrices),sizeof(unsigned int),1,fp);

  /* write attributes */
  fwrite((void *)&(locus->attributes),sizeof(unsigned int),1,fp);

  /* write pattern weights sum */
  fwrite((void *)&(locus->pattern_weights_sum),sizeof(unsigned int),1,fp);

  /* write diploid */
  fwrite((void *)&(locus->diploid),sizeof(int),1,fp);

  if (locus->diploid)
  {
    /* write original diploid number of sites */
    fwrite((void *)&(locus->unphased_length),sizeof(int),1,fp);

    /* write diploid mapping A1 -> A3 */
    fwrite((void *)(locus->diploid_mapping),sizeof(unsigned int),locus->unphased_length,fp);

    /* write diploid resolution count */
    fwrite((void *)(locus->diploid_resolution_count),sizeof(unsigned int),locus->unphased_length,fp);

    /* write pattern weights for original diploid A1 alignment */
    fwrite((void *)(locus->pattern_weights),sizeof(unsigned int),locus->unphased_length,fp);
  }
  else
  {
    fwrite((void *)(locus->pattern_weights),sizeof(unsigned int),locus->sites,fp);
  }

  /* write tip CLVs */
  for (i = 0; i < locus->tips; ++i)
  {
    unsigned int clv_index = gtree->nodes[i]->clv_index;
    long span = locus->sites * locus->states * locus->rate_cats;
    
    fwrite((void *)(locus->clv[clv_index]),sizeof(double),span,fp);
  }
}

static void dump_chk_section_3(FILE * fp, gtree_t ** gtree_list, long msa_count)
{
  long i;

  for (i = 0; i < msa_count; ++i)
  {
    dump_gene_tree(fp,gtree_list[i]);
  }
}

static void dump_chk_section_4(FILE * fp, gtree_t ** gtree_list, locus_t ** locus_list, long msa_count)
{
  long i;

  for (i = 0; i < msa_count; ++i)
  {
    dump_locus(fp,gtree_list[i], locus_list[i]);
  }

}

int checkpoint_dump(stree_t * stree,
                    gtree_t ** gtree_list,
                    locus_t ** locus_list,
                    double * pjump,
                    unsigned long curstep,
                    long ft_round,
                    long mcmc_offset)
{
  FILE * fp;
  char * s = NULL;

  xasprintf(&s, "%s.%ld.chk", opt_outfile, ++opt_checkpoint_current);

  fprintf(stdout,"\n\nWriting checkpoint file %s\n\n",s);

  fp = fopen(s,"w");
  free(s);
  if (!fp)
  {
    fprintf(stderr, "Cannot open file %s for checkpointing...",s);
    return 0;
  }


  /* write checkpoint header */
  dump_chk_header(fp,stree);

  /* write section 1 */
  dump_chk_section_1(fp,stree,pjump,curstep,ft_round,mcmc_offset);

  /* write section 2 */
  dump_chk_section_2(fp,stree);

  /* write section 3 */
  dump_chk_section_3(fp,gtree_list,stree->locus_count);

  /* write section 4 */
  dump_chk_section_4(fp,gtree_list,locus_list,stree->locus_count);

  fclose(fp);
  
  return 1;
}
