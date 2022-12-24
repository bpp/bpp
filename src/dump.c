/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp.h"

#define PROP_COUNT 5
#define GTR_PROP_COUNT 3
#define CLOCK_PROP_COUNT 5

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
  DUMP(&opt_threads,1,fp);
  DUMP(&opt_threads_start,1,fp);
  DUMP(&opt_threads_step,1,fp);
  unsigned int * rng = get_legacy_rndu_array();
  DUMP(rng,opt_threads,fp);

  /* number of sections */
  unsigned int sections = 3;
  DUMP(&sections,1,fp);

  /* compute length of section 1 */
  unsigned long size_section = 0;
  size_section += sizeof(unsigned long);              /* seed */
  size_section += strlen(opt_msafile)+1;              /* msa filename */
  if (opt_mapfile)
    size_section += strlen(opt_mapfile)+1;            /* imap filename */
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
  size_section += sizeof(unsigned long);              /* output file offset */
  if (opt_print_genetrees)
     size_section += opt_locus_count*sizeof(long);    /* gtree file offsets */
  if (opt_print_rates && opt_clock != BPP_CLOCK_GLOBAL)
     size_section += opt_locus_count*sizeof(long);
    

  size_t pjump_size = PROP_COUNT + 1+1 + GTR_PROP_COUNT + CLOCK_PROP_COUNT + !!opt_migration;

  size_section += pjump_size*sizeof(double);          /* pjump */
  size_section += sizeof(long);                       /* dparam_count */
  size_section += sizeof(long);                       /* ft_round_rj*/
  size_section += sizeof(double);                     /* pjump_rj*/
  size_section += sizeof(long);                       /* ft_round_spr */
  size_section += sizeof(long);                       /* ft_round_snl */
  size_section += sizeof(long);                       /* pjump_spr */
  size_section += sizeof(long);                       /* pjump_snl */
  size_section += sizeof(double);                     /* mean_logl */
  size_section += sizeof(double);                     /* mean_root_age */
  size_section += sizeof(double);                     /* mean_root_theta */
  size_section += sizeof(long);                       /* opt_est_locusrate */
  size_section += sizeof(long);                       /* opt_est_heredity */
  size_section += sizeof(double);                     /* opt_heredity_alpha */
  size_section += sizeof(double);                     /* opt_heredity_beta */
  size_section += 4*sizeof(long);                     /* opt_print_* */
  if (opt_est_locusrate || opt_est_heredity == HEREDITY_ESTIMATE)
    size_section += sizeof(double);                   /* locusrate finetune */

  /* TODO: Check which parameters are written for notheta option */


  /* write section 1 size */
  DUMP(&size_section,1,fp);
}

static void dump_chk_section_1(FILE * fp,
                               stree_t * stree,
                               double * pjump,
                               long curstep,
                               long ft_round,
                               long ndspecies,
                               long mcmc_offset,
                               long out_offset,
                               long * gtree_offset,
                               long * rates_offset,
                               long dparam_count,
                               double * posterior,
                               double * pspecies,
                               long dmodels_count,
                               long ft_round_rj,
                               double pjump_rj,
                               long ft_round_spr,
                               long ft_round_snl,
                               long pjump_spr,
                               long pjump_snl,
                               double mean_logl,
                               long * mean_mrate_row,
                               long * mean_mrate_col,
                               long * mean_mrate_round,
                               double * mean_mrate,
                               double * mean_tau,
                               double * mean_theta,
                               double * mean_phi,
                               long mean_mrate_count,
                               long mean_tau_count,
                               long mean_theta_count,
                               long mean_phi_count,
                               int prec_logpg,
                               int prec_logl)
{
  size_t i;
  unsigned int hoffset = stree->tip_count+stree->inner_count;
  unsigned int total_nodes = hoffset+stree->hybrid_count;

  /* write seed */
  DUMP(&opt_seed,1,fp);

  /* write control file */
  DUMP(opt_cfile,strlen(opt_cfile)+1,fp);

  /* write seqfile */
  DUMP(opt_msafile,strlen(opt_msafile)+1,fp);

  /* write constraintfile */
  DUMP(&opt_constraint_count,1,fp);
  if (opt_constraint_count)
    DUMP(opt_constraintfile,strlen(opt_constraintfile)+1,fp);

  /* write whether we will write a map file */
  long mapfile_present = (opt_mapfile ? 1 : 0);
  DUMP(&mapfile_present,1,fp);

  /* write imap file */
  if (mapfile_present)
    DUMP(opt_mapfile,strlen(opt_mapfile)+1,fp);

  /* write outfile */
  DUMP(opt_outfile,strlen(opt_outfile)+1,fp);

  /* write mcmcfile */
  DUMP(opt_mcmcfile,strlen(opt_mcmcfile)+1,fp);

  /* write checkpint info */
  DUMP(&opt_checkpoint,1,fp);
  DUMP(&opt_checkpoint_current,1,fp);
  DUMP(&opt_checkpoint_initial,1,fp);
  DUMP(&opt_checkpoint_step,1,fp);

  /* write network info */
  DUMP(&opt_msci,1,fp);

  /* write migration info */
  DUMP(&opt_migration,1,fp);

  /* write method info */
  DUMP(&opt_method,1,fp);

  /* write speciesdelimitation */
  DUMP(&opt_est_delimit,1,fp);
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
  DUMP(&opt_est_stree,1,fp);
  DUMP(dummy,3*sizeof(double),fp);                   /* typecast */
  
  /* write speciesmodelprior */
  DUMP(&opt_delimit_prior,1,fp);

  /* write species&tree */
  DUMP(&(stree->tip_count),1,fp);
  DUMP(&(stree->inner_count),1,fp);
  DUMP(&(stree->hybrid_count),1,fp);
  DUMP(&(stree->edge_count),1,fp);
  for (i = 0; i < stree->tip_count; ++i)
    DUMP(stree->nodes[i]->label,strlen(stree->nodes[i]->label)+1,fp);
  for (i = 0; i < stree->hybrid_count; ++i)
    DUMP(stree->nodes[hoffset+i]->label,
         strlen(stree->nodes[hoffset+i]->label)+1,
         fp);
  
  /* write usedata, cleandata and nloci */
  DUMP(&opt_usedata,1,fp);
  DUMP(&opt_cleandata,1,fp);
  DUMP(&opt_locus_count,1,fp);

  /* write print flags */
  DUMP(&opt_print_samples,1,fp);
  DUMP(&opt_print_locusrate,1,fp);
  DUMP(&opt_print_hscalars,1,fp);
  DUMP(&opt_print_genetrees,1,fp);
  DUMP(&opt_print_rates,1,fp);
  DUMP(&opt_print_qmatrix,1,fp);
  DUMP(&opt_print_locusfile,1,fp);

  /* write theta prior */
  DUMP(&opt_theta_dist,1,fp);
  DUMP(&opt_theta_alpha,1,fp);
  DUMP(&opt_theta_beta,1,fp);
  DUMP(&opt_theta_p,1,fp);
  DUMP(&opt_theta_q,1,fp);
  DUMP(&opt_theta_min,1,fp);
  DUMP(&opt_theta_max,1,fp);
  DUMP(&opt_est_theta,1,fp);
  DUMP(&opt_linkedtheta,1,fp);

  /* write tau prior */
  DUMP(&opt_tau_dist,1,fp);
  DUMP(&opt_tau_alpha,1,fp);
  DUMP(&opt_tau_beta,1,fp);

  DUMP(&opt_phi_alpha,1,fp);
  DUMP(&opt_phi_beta,1,fp);

  DUMP(&opt_mig_alpha,1,fp);
  DUMP(&opt_mig_beta,1,fp);

  /* write substitution model information */
  DUMP(&opt_model,1,fp);

  /* write gamma rate variation information */
  DUMP(&opt_alpha_cats,1,fp);
  DUMP(&opt_alpha_alpha,1,fp);
  DUMP(&opt_alpha_beta,1,fp);


  /* whether locus mutation rate is estimated */
  DUMP(&opt_est_locusrate,1,fp);

  /* whether mubar is estimated */
  DUMP(&opt_est_mubar,1,fp);

  /* whether heredity scalers are estimated */
  DUMP(&opt_est_heredity,1,fp);
  DUMP(&opt_heredity_alpha,1,fp);
  DUMP(&opt_heredity_beta,1,fp);

  /* write clock and locusrate info */
  DUMP(&opt_clock,1,fp);
  DUMP(&opt_mubar_alpha,1,fp);
  DUMP(&opt_mubar_beta,1,fp);
  DUMP(&opt_mui_alpha,1,fp);
  DUMP(&opt_vbar_alpha,1,fp);
  DUMP(&opt_vbar_beta,1,fp);
  DUMP(&opt_vi_alpha,1,fp);
  DUMP(&opt_rate_prior,1,fp);
  DUMP(&opt_locusrate_prior,1,fp);

  /* write finetune */
  DUMP(&opt_finetune_reset,1,fp);
  DUMP(&opt_finetune_migrates,1,fp);
  DUMP(&opt_finetune_phi,1,fp);
  DUMP(&opt_finetune_gtage,1,fp);
  DUMP(&opt_finetune_gtspr,1,fp);
  DUMP(&opt_finetune_theta,1,fp);
  DUMP(&opt_finetune_tau,1,fp);
  DUMP(&opt_finetune_mix,1,fp);
  DUMP(&opt_finetune_locusrate,1,fp);
  DUMP(&opt_finetune_qrates,1,fp);
  DUMP(&opt_finetune_freqs,1,fp);
  DUMP(&opt_finetune_alpha,1,fp);
  DUMP(&opt_finetune_mubar,1,fp);
  DUMP(&opt_finetune_mui,1,fp);
  DUMP(&opt_finetune_nubar,1,fp);
  DUMP(&opt_finetune_nui,1,fp);
  DUMP(&opt_finetune_branchrate,1,fp);

  DUMP(&opt_max_species_count,1,fp);

  DUMP(&opt_prob_snl,1,fp);
  DUMP(&opt_prob_snl_shrink,1,fp);
  DUMP(&opt_snl_lambda_expand,1,fp);
  DUMP(&opt_snl_lambda_shrink,1,fp);

  DUMP(&opt_pseudop_exist,1,fp);
  DUMP(&opt_mig_vrates_exist,1,fp);

  DUMP(&(stree->locusrate_mubar),1,fp);
  DUMP(&(stree->locusrate_nubar),1,fp);
  DUMP(&(stree->nui_sum),1,fp);

  /* write diploid */
  if (opt_diploid)
    DUMP(opt_diploid,stree->tip_count,fp);
  else
  {
    for (i = 0; i < stree->tip_count; ++i)
      DUMP(dummy,sizeof(long),fp);         /* typecast */
  }
  DUMP(&opt_diploid_size,1,fp);

  /* write mcmc run info */
  DUMP(&opt_burnin,1,fp);
  DUMP(&opt_samplefreq,1,fp);
  DUMP(&opt_samples,1,fp);
  DUMP(&curstep,1,fp);
  DUMP(&ft_round,1,fp);
  DUMP(&ndspecies,1,fp);

  size_t pjump_size = PROP_COUNT + 1+1 + GTR_PROP_COUNT + CLOCK_PROP_COUNT + !!opt_migration;
  /* write pjump */
  DUMP(pjump,pjump_size,fp);

  /* write MCMC file offset */
  DUMP(&mcmc_offset,1,fp);

  /* write output file offset */
  DUMP(&out_offset,1,fp);

  /* write bfbeta */
  DUMP(&opt_bfbeta,1,fp);

  /* write gtree file offset if available*/
  if (opt_print_genetrees)
    DUMP(gtree_offset,opt_locus_count,fp);

  if (opt_print_locusfile)
    DUMP(rates_offset,opt_locus_count,fp);

  DUMP(&dparam_count,1,fp);

  DUMP(&dmodels_count,1,fp);
  if (dmodels_count)
    DUMP(posterior,dmodels_count,fp);

  if (opt_method == METHOD_11)
    DUMP(pspecies,opt_max_species_count,fp);

  DUMP(&ft_round_rj,1,fp);
  DUMP(&pjump_rj,1,fp);
  DUMP(&ft_round_spr,1,fp);
  DUMP(&ft_round_snl, 1, fp);
  DUMP(&pjump_spr, 1, fp);
  DUMP(&pjump_snl,1,fp);
  DUMP(&mean_logl,1,fp);
  if (opt_migration)
  {
    DUMP(&mean_mrate_count,1,fp);
    DUMP(mean_mrate,mean_mrate_count,fp);
    DUMP(mean_mrate_row, mean_mrate_count, fp);
    DUMP(mean_mrate_col, mean_mrate_count, fp);
    DUMP(mean_mrate_round, mean_mrate_count, fp);
  }
  DUMP(&mean_tau_count,1,fp);
  if (opt_est_theta)
    DUMP(&mean_theta_count,1,fp);
  DUMP(&mean_phi_count,1,fp);

  DUMP(mean_tau,mean_tau_count,fp);
  if (opt_est_theta)
    DUMP(mean_theta,mean_theta_count,fp);

  if (mean_phi_count)
  {
    assert(opt_msci);
    DUMP(mean_phi,mean_phi_count,fp);
  }

  DUMP(&prec_logpg,1,fp);
  DUMP(&prec_logl,1,fp);

  DUMP(&opt_load_balance,1,fp);

  if (opt_threads > 1)
  {
    thread_info_t * ti = threads_ti();
    for (i = 0; i < (size_t)opt_threads; ++i)
    {
      thread_info_t * tip = ti+i;
      DUMP(&(tip->locus_first),1,fp);
      DUMP(&(tip->locus_count),1,fp);
    }
  }

  if (opt_migration)
  {
    assert(stree->hybrid_count == 0);
    for (i = 0; i < total_nodes; ++i)
      DUMP(opt_migration_matrix[i],total_nodes,fp);

    for (i = 0; i < total_nodes; ++i)
      DUMP(opt_mig_bitmatrix[i],total_nodes,fp);

    for (i = 0; i < opt_migration; ++i)
    {
      migspec_t * spec = opt_mig_specs+i;
      DUMP(spec->source, strlen(spec->source)+1, fp);
      DUMP(spec->target, strlen(spec->target)+1, fp);
      DUMP(&(spec->si), 1, fp);
      DUMP(&(spec->ti), 1, fp);
      DUMP(&(spec->am), 1, fp);
      DUMP(&(spec->alpha), 1, fp);
      DUMP(&(spec->beta), 1, fp);
      DUMP(&(spec->pseudo_a), 1, fp);
      DUMP(&(spec->pseudo_b), 1, fp);
      DUMP(&(spec->params), 1, fp);
      DUMP(&(spec->M), 1, fp);
      if (spec->params == 1 || spec->params == 3 || spec->params == 5)
        DUMP(spec->Mi, opt_locus_count, fp);
    }
  }
}


static void dump_chk_section_2(FILE * fp, stree_t * stree)
{
  unsigned int total_nodes;
  unsigned int hoffset;
  long i,j;
  long ltheta_index;

  total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  /* write node indices of hybridization events */
  hoffset = stree->tip_count + stree->inner_count;
  for (i = 0; i < stree->hybrid_count; ++i)
  {
    assert(node_is_mirror(stree->nodes[hoffset+i]));
    DUMP(&(stree->nodes[hoffset+i]->hybrid->node_index),1,fp);
  }

  /* write left child node indices */
  for (i = 0; i < stree->inner_count; ++i)
    DUMP(&(stree->nodes[stree->tip_count+i]->left->node_index),1,fp);


  /* write right child node indices */
  unsigned int valid;
  for (i = 0; i < stree->inner_count; ++i)
  {
    if (stree->nodes[stree->tip_count+i]->right)
    {
      valid = 1;
      DUMP(&valid,1,fp);
      DUMP(&(stree->nodes[stree->tip_count+i]->right->node_index),1,fp);
    }
    else
    {
      valid = 0;
      DUMP(&valid,1,fp);
    }
  }

  for (i = 0; i < stree->hybrid_count; ++i)
    DUMP(&(stree->nodes[hoffset+i]->hybrid->hphi),1,fp);

  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->has_phi),1,fp);

  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->htau),1,fp);

  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->prop_tau),1,fp);

  for (i = 0; i < total_nodes; ++i)
  {
    ltheta_index = stree->nodes[i]->linked_theta ?
                     (long)stree->nodes[i]->linked_theta->node_index : -1;
    DUMP(&ltheta_index,1,fp);
  }


  /* TODO: We do not need to write theta when !opt_est_theta */
  /* write theta */
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->theta),1,fp);
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->has_theta),1,fp);

  /* write tau */
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->tau),1,fp);

  /* write support */
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->support),1,fp);

  /* write constraints */
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->constraint),1,fp);

  /* write constraints line numbers */
  for (i = 0; i < total_nodes; ++i)
    DUMP(&(stree->nodes[i]->constraint_lineno),1,fp);

  /* write number of coalescent events */
  assert(opt_locus_count == stree->locus_count);
  for (i = 0; i < total_nodes; ++i)
    DUMP(stree->nodes[i]->event_count,opt_locus_count,fp);

  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    for (i = 0; i < total_nodes; ++i)
    {
      if (stree->nodes[i]->brate)
      {
        valid = 1;
        DUMP(&valid,1,fp);
        DUMP(stree->nodes[i]->brate,opt_locus_count,fp);
      }
      else
      {
        valid = 0;
        DUMP(&valid,1,fp);
      }
    }
  }

//  /* write MSC density contribution - TODO: Candidate for removal */
//  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
//    fwrite((void *)&(stree->nodes[i]->logpr_contrib),sizeof(double),1,fp);
  
  /* TODO: Perhaps we can remove this and compute from scratch when resuming */
  if (!opt_est_theta)
  {
    DUMP(&(stree->notheta_logpr),1,fp);
    DUMP(&(stree->notheta_hfactor),1,fp);
    DUMP(&(stree->notheta_sfactor),1,fp);
    for (i = 0; i < total_nodes; ++i)
    {
      DUMP(stree->nodes[i]->t2h,opt_locus_count,fp);
      DUMP(&(stree->nodes[i]->t2h_sum),1,fp);
      DUMP(&(stree->nodes[i]->event_count_sum),1,fp);
      DUMP(&(stree->nodes[i]->notheta_logpr_contrib),1,fp);
    }

    if (opt_msci)
    {
      for (i = 0; i < stree->hybrid_count; ++i)
      {
        unsigned int index = stree->tip_count+stree->inner_count;
        snode_t * x = stree->nodes[index+i];

        DUMP(x->notheta_phi_contrib, opt_locus_count, fp);
        DUMP(x->hybrid->notheta_phi_contrib, opt_locus_count, fp);
        DUMP(&(x->hphi_sum),1,fp);
        DUMP(&(x->hybrid->hphi_sum),1,fp);
      }
    }
  }

  DUMP(&(stree->root_age),1,fp);

  /* TODO : Perhaps write only seqin_count for tips? */
  /* write number of incoming sequences for each node */
  for (i = 0; i < total_nodes; ++i)
    DUMP(stree->nodes[i]->seqin_count,opt_locus_count,fp);

  /* write event indices for each node */
  for (i = 0; i < total_nodes; ++i)
  {
    for (j = 0; j < opt_locus_count; ++j)
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

  /* write migevent_count */
  if (opt_migration)
  {
    for (i = 0; i < total_nodes; ++i)
      DUMP(stree->nodes[i]->migevent_count,opt_locus_count,fp);

    for (i = 0; i < total_nodes; ++i)
    {
      DUMP(&(stree->nodes[i]->mb_count),1,fp);
      DUMP(stree->nodes[i]->migbuffer,stree->nodes[i]->mb_count,fp);
    }
  }
}

static void dump_gene_tree(FILE * fp, gtree_t * gtree, stree_t * stree)
{
  long i,j;

  unsigned int hybrid_count = stree->hybrid_count;

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

  /* write hpath */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    DUMP(gtree->nodes[i]->hpath,hybrid_count,fp);

  DUMP(&(gtree->rate_mui),1,fp);
  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    DUMP(&(gtree->rate_nui),1,fp);
    DUMP(&(gtree->lnprior_rates),1,fp);
  }

  DUMP(&(gtree->original_index),1,fp);

  DUMP(&(gtree->msa_index),1,fp);

  if (opt_migration)
  {
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    {
      miginfo_t * mi = gtree->nodes[i]->mi;

      if (!mi || mi->count == 0)
      {
        long mi_count = 0;
        DUMP(&mi_count,1,fp);
        continue;
      }

      DUMP(&(mi->count),1,fp);
      for (j = 0; j < mi->count; ++j)
      {
        DUMP(&(mi->me[j].time),1,fp);
        DUMP(&(mi->me[j].source->node_index),1,fp);
        DUMP(&(mi->me[j].target->node_index),1,fp);
      }
    }

    /* migcount */
    long total_snodes = stree->tip_count+stree->inner_count;
    for (i = 0; i < total_snodes; ++i)
      DUMP(gtree->migcount[i],total_snodes,fp);
  }
  
}

static void dump_locus(FILE * fp, gtree_t * gtree, locus_t * locus)
{

  long i;

  /* write data type */
  DUMP(&(locus->dtype),1,fp);

  /* write substitution model */
  DUMP(&(locus->model),1,fp);

  /* write number of sites */
  DUMP(&(locus->sites),1,fp);

  /* write number of states */
  DUMP(&(locus->states),1,fp);

  /* write number of rate categories */
  DUMP(&(locus->rate_cats),1,fp);

  /* write number of rate matrices */
  DUMP(&(locus->rate_matrices),1,fp);

  /* write number of prob matrices */
  DUMP(&(locus->prob_matrices),1,fp);

  /* write number of prob matrices */
  DUMP(&(locus->scale_buffers),1,fp);

  /* write attributes */
  DUMP(&(locus->attributes),1,fp);

  /* write pattern weights sum */
  DUMP(&(locus->pattern_weights_sum),1,fp);

  /* write alpha */
  DUMP(&(locus->rates_alpha),1,fp);

  /* write qrates param count */
  DUMP(&(locus->qrates_param_count),1,fp);

  /* write freqs param count */
  DUMP(&(locus->freqs_param_count),1,fp);

  /* write category rates */
  DUMP(locus->rates,locus->rate_cats,fp);

  /* dump base frequencies */
  for (i = 0; i < locus->rate_matrices; ++i)
    DUMP(locus->frequencies[i],locus->states,fp);

  /* dump qmatrix rates */
  for (i = 0; i < locus->rate_matrices; ++i)
    DUMP(locus->subst_params[i],((locus->states-1)*locus->states)/2,fp);

  /* write param indices */
  DUMP(locus->param_indices,locus->rate_cats,fp);

  /* write heredity scalars */
  DUMP(locus->heredity,locus->rate_matrices,fp);

  /* write diploid */
  DUMP(&(locus->diploid),1,fp);

  if (locus->diploid)
  {
    size_t sites_a2 = 0;

    /* write original diploid number of sites */
    DUMP(&(locus->unphased_length),1,fp);

    /* write diploid resolution count (A1 -> A2)*/
    DUMP(locus->diploid_resolution_count,locus->unphased_length,fp);

    for (i = 0; i < locus->unphased_length; ++i)
      sites_a2 += locus->diploid_resolution_count[i];

    /* write diploid mapping A2 -> A3 */
    DUMP(locus->diploid_mapping,sites_a2,fp);

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

  DUMP(&(locus->original_index),1,fp);
}

static void dump_chk_section_3(FILE * fp, gtree_t ** gtree_list, stree_t * stree, long msa_count)
{
  long i;

  for (i = 0; i < msa_count; ++i)
  {
    dump_gene_tree(fp,gtree_list[i],stree);
  }
}

static void dump_chk_section_4(FILE * fp,
                               gtree_t ** gtree_list,
                               locus_t ** locus_list,
                               long msa_count)
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
                    long ndspecies,
                    long mcmc_offset,
                    long out_offset,
                    long * gtree_offset,
                    long * rates_offset,
                    long dparam_count,
                    double * posterior,
                    double * pspecies,
                    long dmodels_count,
                    long ft_round_rj,
                    double pjump_rj,
                    long ft_round_spr,
                    long ft_round_snl,
                    long pjump_spr,
                    long pjump_snl,
                    double mean_logl,
                    long * mean_mrate_row,
                    long * mean_mrate_col,
                    long * mean_mrate_round,
                    double * mean_mrate,
                    double * mean_tau,
                    double * mean_theta,
                    double * mean_phi,
                    long mean_mrate_count,
                    long mean_tau_count,
                    long mean_theta_count,
                    long mean_phi_count,
                    int prec_logpg,
                    int prec_logl)
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
  dump_chk_section_1(fp,
                     stree,
                     pjump,
                     curstep,
                     ft_round,
                     ndspecies,
                     mcmc_offset,
                     out_offset,
                     gtree_offset,
                     rates_offset,
                     dparam_count,
                     posterior,
                     pspecies,
                     dmodels_count,
                     ft_round_rj,
                     pjump_rj,
                     ft_round_spr,
                     ft_round_snl,
                     pjump_spr,
                     pjump_snl,
                     mean_logl,
                     mean_mrate_row,
                     mean_mrate_col,
                     mean_mrate_round,
                     mean_mrate,
                     mean_tau,
                     mean_theta,
                     mean_phi,
                     mean_mrate_count,
                     mean_tau_count,
                     mean_theta_count,
                     mean_phi_count,
                     prec_logpg,
                     prec_logl);

  /* write section 2 */
  dump_chk_section_2(fp,stree);

  /* write section 3 */
  dump_chk_section_3(fp,gtree_list,stree,stree->locus_count);

  /* write section 4 */
  dump_chk_section_4(fp,gtree_list,locus_list,stree->locus_count);

  fclose(fp);
  
  return 1;
}
