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

#define DUMP(x,n,fp) fwrite((void *)(x),sizeof(*(x)),n,fp)

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

  DUMP(header,20,fp);

  size_t size_int = sizeof(int);
  size_t size_long = sizeof(long);
  size_t size_double = sizeof(double);

  /* write size of int */
  size_type = (BYTE)(size_int & 0xFF);
  DUMP(&size_type,1,fp);

  /* write size of long */
  size_type = (BYTE)(size_long & 0xFF);
  DUMP(&size_type,1,fp);

  /* write size of double */
  size_type = (BYTE)(size_double & 0xFF);
  DUMP(&size_type,1,fp);

  /* write RNG value */
  unsigned int rng = get_legacy_rndu_status();
  DUMP(&rng,1,fp);

  /* number of sections */
  unsigned int sections = 3;
  DUMP(&sections,1,fp);

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
  DUMP(&size_section,1,fp);
}

static void dump_chk_section_1(FILE * fp, stree_t * stree, double * pjump, long curstep, long ft_round, long mcmc_offset)
{
  size_t i;

  /* write seed */
  DUMP(&opt_seed,1,fp);

  /* write seqfile */
  DUMP(opt_msafile,strlen(opt_msafile)+1,fp);

  /* write imap file */
  DUMP(opt_mapfile,strlen(opt_mapfile)+1,fp);

  /* write outfile */
  DUMP(opt_outfile,strlen(opt_outfile)+1,fp);

  /* write mcmcfile */
  DUMP(opt_mcmcfile,strlen(opt_mcmcfile)+1,fp);

  /* write speciesdelimitation */
  DUMP(&opt_delimit,1,fp);
  DUMP(&opt_rjmcmc_method,1,fp);

  /* always write two more elements, even if speciesdelimitation = 0 */
  if (opt_rjmcmc_method == 0)
  {
    DUMP(&opt_rjmcmc_epsilon,1,fp);

    /* dummy bytes to fill the file */
    DUMP(dummy,sizeof(double),fp);              /* typecast */
  }
  else
  {
    DUMP(&opt_rjmcmc_alpha,1,fp);
    DUMP(&opt_rjmcmc_mean,1,fp);
  }

  /* write speciestree */
  DUMP(&opt_stree,1,fp);
  DUMP(dummy,3*sizeof(double),fp);                   /* typecast */
  
  /* write speciesmodelprior */
  DUMP(&opt_delimit_prior,1,fp);

  /* write species&tree */
  DUMP(&(stree->tip_count),1,fp);
  for (i = 0; i < stree->tip_count; ++i)
    DUMP(stree->nodes[i]->label,strlen(stree->nodes[i]->label)+1,fp);
  
  /* write usedata, cleandata and nloci */
  DUMP(&opt_usedata,1,fp);
  DUMP(&opt_cleandata,1,fp);
  DUMP(&opt_nloci,1,fp);

  /* write theta prior */
  DUMP(&opt_theta_alpha,1,fp);
  DUMP(&opt_theta_beta,1,fp);
  DUMP(&opt_est_theta,1,fp);

  /* write tau prior */
  DUMP(&opt_tau_alpha,1,fp);
  DUMP(&opt_tau_beta,1,fp);

  /* write finetune */
  DUMP(&opt_finetune_reset,1,fp);
  DUMP(&opt_finetune_gtage,1,fp);
  DUMP(&opt_finetune_gtspr,1,fp);
  DUMP(&opt_finetune_theta,1,fp);
  DUMP(&opt_finetune_tau,1,fp);
  DUMP(&opt_finetune_mix,1,fp);

  /* write diploid */
  if (opt_diploid)
    DUMP(opt_diploid,stree->tip_count,fp);
  else
  {
    for (i = 0; i < stree->tip_count; ++i)
      DUMP(dummy,sizeof(long),fp);         /* typecast */
  }

  /* write mcmc run info */
  DUMP(&opt_burnin,1,fp);
  DUMP(&opt_samplefreq,1,fp);
  DUMP(&opt_samples,1,fp);
  DUMP(&curstep,1,fp);
  DUMP(&ft_round,1,fp);

  /* write pjump */
  DUMP(pjump,PROP_COUNT,fp);

  /* write MCMC file offset */
  DUMP(&mcmc_offset,1,fp);
}


static void dump_chk_section_2(FILE * fp, stree_t * stree)
{
  long i,j;

  /* write left child node indices */
  for (i = 0; i < stree->inner_count; ++i)
    DUMP(&(stree->nodes[stree->tip_count+i]->left->node_index),1,fp);

  /* write right child node indices */
  for (i = 0; i < stree->inner_count; ++i)
    DUMP(&(stree->nodes[stree->tip_count+i]->right->node_index),1,fp);

  /* write theta */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    DUMP(&(stree->nodes[i]->theta),1,fp);

  /* write tau */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    DUMP(&(stree->nodes[i]->tau),1,fp);

  /* write support */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    DUMP(&(stree->nodes[i]->support),1,fp);

  /* write number of coalescent events */
  assert(opt_nloci == stree->locus_count);
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    DUMP(stree->nodes[i]->event_count,opt_nloci,fp);

//  /* write MSC density contribution - TODO: Candidate for removal */
//  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
//    fwrite((void *)&(stree->nodes[i]->logpr_contrib),sizeof(double),1,fp);

  DUMP(&(stree->root_age),1,fp);

  /* TODO : Perhaps write only seqin_count for tips? */
  /* write number of incoming sequences for each node */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    DUMP(stree->nodes[i]->seqin_count,opt_nloci,fp);

  /* write event indices for each node */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    for (j = 0; j < opt_nloci; ++j)
    {
      dlist_item_t * di = stree->nodes[i]->event[j]->head;
      while (di)
      {
        gnode_t * gnode = (gnode_t *)(di->data);

        DUMP(&(gnode->node_index),1,fp);

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
    DUMP(gtree->nodes[i]->label,strlen(gtree->nodes[i]->label)+1,fp);

  /* write left child node indices */
  for (i = 0; i < gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[gtree->tip_count+i]->left->node_index),1,fp);

  /* write right child node indices */
  for (i = 0; i < gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[gtree->tip_count+i]->right->node_index),1,fp);

  /* write branch lengths - TODO: Candidate for removal */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->length),1,fp);

  /* write ages */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->time),1,fp);

  /* write population index (corresponding species tree node index) */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->pop->node_index),1,fp);

  /* write CLV indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->clv_index),1,fp);

  /* write scaler indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->scaler_index),1,fp);

  /* write pmatrix indices */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->pmatrix_index),1,fp);

  /* write mark - TODO: Candidate for removal */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(&(gtree->nodes[i]->mark),1,fp);
}

static void dump_locus(FILE * fp, gtree_t * gtree, locus_t * locus)
{

  long i;

  /* write number of sites */
  DUMP(&(locus->sites),1,fp);

  /* write number of states */
  DUMP(&(locus->states),1,fp);

  /* write number of rate categories */
  DUMP(&(locus->rate_cats),1,fp);

  /* write number of rate matrices */
  DUMP(&(locus->rate_matrices),1,fp);

  /* write attributes */
  DUMP(&(locus->attributes),1,fp);

  /* write pattern weights sum */
  DUMP(&(locus->pattern_weights_sum),1,fp);

  /* write diploid */
  DUMP(&(locus->diploid),1,fp);

  if (locus->diploid)
  {
    /* write original diploid number of sites */
    DUMP(&(locus->unphased_length),1,fp);

    /* write diploid mapping A1 -> A3 */
    DUMP(locus->diploid_mapping,locus->unphased_length,fp);

    /* write diploid resolution count */
    DUMP(locus->diploid_resolution_count,locus->unphased_length,fp);

    /* write pattern weights for original diploid A1 alignment */
    DUMP(locus->pattern_weights,locus->unphased_length,fp);
  }
  else
  {
    DUMP(locus->pattern_weights,locus->sites,fp);
  }

  /* write tip CLVs */
  for (i = 0; i < locus->tips; ++i)
  {
    unsigned int clv_index = gtree->nodes[i]->clv_index;
    long span = locus->sites * locus->states * locus->rate_cats;
    
    DUMP(locus->clv[clv_index],span,fp);
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