/*
    Copyright (C) 2016-2017 Tomas Flouri and Ziheng Yang

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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))

#define MARK_BRANCH_UPDATE      1
#define MARK_ANCESTOR_LNODE     2
#define MARK_ANCESTOR_RNODE     4
#define MARK_POP_CHANGE         8
#define MARK_AGE_UPDATE        16

#define lbeta(p,q) (lgamma(p) + lgamma(q) - lgamma(p+q))
#define log_pdfbeta(x,p,q,b) (-lbeta(p,q) + (p-1)*log(x/b) + (q-1)*log(1-x/b) -\
        log(b))
#define log_pdfinvgamma(x, a, b)  ( (a)*log(b) - lgamma(a) - (a+1)*log(x) - (b)/(x) )

static int opt_rjalgorithm = 1;
static double opt_finetune_rja = 2;
static double opt_finetune_rjm = 1;

static gnode_t ** nodevec;
static unsigned int * nodevec_offset;
static unsigned int * nodevec_count;

static int * feasible;
static gnode_t *** partials;
static unsigned int * partials_count;

static int dbg_rj_join = 0;

#if 0
static void all_partials_recursive(gnode_t * node,
                                   unsigned int * trav_size,
                                   gnode_t ** outbuffer)
{
  if (!node->left)
    return;

  all_partials_recursive(node->left,  trav_size, outbuffer);
  all_partials_recursive(node->right, trav_size, outbuffer);

  outbuffer[*trav_size] = node;
  *trav_size = *trav_size + 1;
}

static void gtree_all_partials(gnode_t * root,
                               gnode_t ** travbuffer,
                               unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return;

  all_partials_recursive(root, trav_size, travbuffer);
}
#endif

void rj_init(gtree_t ** gtreelist, stree_t * stree, unsigned int count)
{
  unsigned int i;
  unsigned int sum = 0;

  for (i = 0; i < count; ++i)
    sum += gtreelist[i]->tip_count + gtreelist[i]->inner_count;

  nodevec = (gnode_t **)xcalloc(sum,sizeof(gnode_t *));
  nodevec_count = (unsigned int *)xcalloc(count,sizeof(unsigned int));

  /* get offset to the space allocated to each gene tree for storing nodes */
  nodevec_offset = (unsigned int *)xmalloc(count*sizeof(unsigned int));
  for (sum = 0, i = 0; i < count; ++i)
  {
    nodevec_offset[i] = sum;
    sum += gtreelist[i]->tip_count + gtreelist[i]->inner_count;
  }
  feasible = (int *)xmalloc((stree->tip_count + stree->inner_count) *
                            sizeof(int));
  partials = (gnode_t ***)xmalloc(stree->locus_count * sizeof(gnode_t **));
  partials_count = (unsigned int *)xcalloc(stree->locus_count,sizeof(unsigned int));
}
void rj_fini()
{
  free(nodevec);
  free(nodevec_offset);
  free(nodevec_count);
  free(feasible);
}

/* TODO : this is the same as gtree.c */
static void unlink_event(gnode_t * node, int msa_index)
{
  /* first re-link the event before the current node with the one after */
  if (node->event->prev)
    node->event->prev->next = node->event->next;
  else
    node->pop->event[msa_index]->head = node->event->next;

  /* now re-link the event after the current node with the one before */
  if (node->event->next)
    node->event->next->prev = node->event->prev;
  else
    node->pop->event[msa_index]->tail = node->event->prev;
}
static double pdf_gamma(double x, double alpha, double beta)
{
/* gamma density: mean=alpha/beta; var=alpha/beta^2
*/
   if (x<=0 || alpha<=0 || beta<=0) {
      printf("x=%.6f a=%.6f b=%.6f", x, alpha, beta);
      fatal("x a b outside range in PDFGamma()");
   }
   if (alpha>100)
      fatal("large alpha in PDFGamma()");
   return pow(beta*x,alpha)/x * exp(-beta*x - lgamma(alpha));
}

static void locate_nodes(stree_t * stree,
                         snode_t * snode,
                         gtree_t * gtree,
                         double tau_upper,
                         int msa_index)
{
  unsigned int i;
  gnode_t * x;
  gnode_t * gnode;

  /* mark all coalescent events (inner nodes) that belong to population snode
     and have two lineages */
  /* TODO: Sanity check - Remove */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
    assert(gtree->nodes[i]->mark == 0);

  /* mark all nodes that are on a lineage that passes left ancestor, upto tau
     upper, which cannot exceed population snode */
  for (i = 0; i < gtree->tip_count; ++i)
  {
    gnode = gtree->nodes[i];
    if (stree->pptable[gnode->pop->node_index][snode->left->node_index])
    {
      gnode->mark = MARK_ANCESTOR_LNODE;

      /* TODO: Mark should be checked for ANCESTOR_(L/R)NODE only and not 0 */
      for (x=gnode->parent; x && x->mark==0 && x->time<=tau_upper; x=x->parent)
        x->mark = MARK_ANCESTOR_LNODE;   /* TODO : We should OR it */
    }
  }

  /* now mark all nodes that are on a lineage that passes right ancestor, upto
     tau upper, which cannot exceed population snode */
  for (i = 0; i < gtree->tip_count; ++i)
  {
    gnode = gtree->nodes[i];
    if (stree->pptable[gnode->pop->node_index][snode->right->node_index])
    {
      gnode->mark = MARK_ANCESTOR_RNODE;     /* TODO: we should OR it */

      for (x=gnode->parent; x && !(x->mark & MARK_ANCESTOR_RNODE) && x->time<=tau_upper; x=x->parent)
        x->mark |= MARK_ANCESTOR_RNODE;
    }
  }
}

static int rubber_update(gnode_t * node, double term, int msa_index)
{
  int nwithin = 1;
  gnode_t ** nv_locus = nodevec + nodevec_offset[msa_index];

  if (!node->left)
  {
    if ((node->mark & MARK_BRANCH_UPDATE) == 0)
    {
      node->mark = MARK_BRANCH_UPDATE;
      nv_locus[nodevec_count[msa_index]++] = node;
    }
    
    return 0;
  }

  /* This does not hold
  assert((node->mark & MARK_BRANCH_UPDATE) == 0);
  */

  if ((node->mark & MARK_BRANCH_UPDATE) == 0)
  {
    node->mark |= MARK_BRANCH_UPDATE;
    nv_locus[nodevec_count[msa_index]++] = node;
  }

  node->old_time = node->time;
  node->time *= term;
  node->mark |= MARK_AGE_UPDATE;

  nwithin += rubber_update(node->left,term,msa_index);
  nwithin += rubber_update(node->right,term,msa_index);

  return nwithin;


  /* children will be marked in recursion */
}

static int rubber_proportional(stree_t * stree,
                               snode_t * snode,
                               gtree_t ** gtreelist,
                               double tau_upper,
                               double tau,
                               double tau_new,
                               int msa_index,
                               double * lnacceptance)
{
  unsigned int i;
  int changed_count = 0;
  int nwithin;
  double rubber = (tau_upper - tau_new) / (tau_upper - tau);
  double tnew;
  double y;

  gtree_t * gtree = gtreelist[msa_index];
  gnode_t ** nv_locus = nodevec + nodevec_offset[msa_index];

  /* mark nodes */
  locate_nodes(stree,snode,gtree,tau_upper,msa_index);

  /* Go through all nodes of the snode population that have lineages coming
     from both child populations and have time <= tau_upper */
  dlist_item_t * event;
  for (event = snode->event[msa_index]->head; event; event = event->next)
  {
    gnode_t * tmp = (gnode_t *)(event->data);

    if ((tmp->mark & (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE)) != 
        (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE))
      continue;

    changed_count++;

    tnew = tau_upper - rubber*(tau_upper - tmp->time);

    tmp->old_time = tmp->time;

    tmp->time = tnew;
    tmp->mark |= MARK_AGE_UPDATE;

    /* add nodes whose branch lengths will change into the list */
    if (!(tmp->mark & MARK_BRANCH_UPDATE))
    {
      if (tmp->parent)
        nv_locus[nodevec_count[msa_index]++] = tmp;

      tmp->mark |= MARK_BRANCH_UPDATE;
    }

    if (!(tmp->left->mark & MARK_BRANCH_UPDATE))
    {
      nv_locus[nodevec_count[msa_index]++] = tmp->left;

      tmp->left->mark |= MARK_BRANCH_UPDATE;
    }

    if (!(tmp->right->mark & MARK_BRANCH_UPDATE))
    {
      nv_locus[nodevec_count[msa_index]++] = tmp->right;

      tmp->right->mark |= MARK_BRANCH_UPDATE;
    }

  }

  /* if node schanged, then update the age of all descendant nodes */
  nwithin = -1;
  y = -1;
  if (changed_count)
  {
    y = 1;
    nwithin = 0;

    for (event = snode->event[msa_index]->head; event; event = event->next)
    {
      gnode_t * tmp = (gnode_t *)(event->data);

      if ((tmp->mark & (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE)) !=
          (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE))
        continue;

      unsigned int count = 0;

      if ((tmp->left->mark & (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE)) !=
          (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE))
        count += rubber_update(tmp->left, tmp->time/tmp->old_time, msa_index);

      if ((tmp->right->mark & (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE)) !=
          (MARK_ANCESTOR_LNODE | MARK_ANCESTOR_RNODE))
        count += rubber_update(tmp->right, tmp->time/tmp->old_time, msa_index);

      for (i = 0; i < count; ++i)
        y *= (tmp->time/tmp->old_time);

      nwithin += count;
    }

    if (nwithin)
      *lnacceptance += log(y);
    *lnacceptance += changed_count * log(rubber);
  }

  if (tau_new > 0)
  {
    /* this is done for SPLIT */
    for (i = gtree->tip_count; i < gtree->tip_count + gtree->inner_count; ++i)
    {
      gnode_t * gnode = gtree->nodes[i];

      if (gnode->pop == snode && gnode->time < tau_new)
      {
        snode_t * newpop;
        if (gnode->mark & MARK_ANCESTOR_LNODE)
          newpop = snode->left;
        else
          newpop = snode->right;
        
        unlink_event(gnode, msa_index);

        gnode->pop->event_count[msa_index]--;

        gnode->mark |= MARK_POP_CHANGE;

        gnode->pop = newpop;

        dlist_item_append(gnode->pop->event[msa_index], gnode->event);

        gnode->pop->event_count[msa_index]++;
        
        snode->seqin_count[msa_index]--;
      }
    }
  }
  else
  {
    /* this is done for JOIN */
    for (i = gtree->tip_count; i < gtree->tip_count + gtree->inner_count; ++i)
    {
      gnode_t * gnode = gtree->nodes[i];

      if (gnode->pop == snode->left || gnode->pop == snode->right)
      {
        unlink_event(gnode, msa_index);

        gnode->pop->event_count[msa_index]--;

        gnode->mark |= MARK_POP_CHANGE;

        gnode->old_pop = gnode->pop;
        gnode->pop = snode;

        dlist_item_append(gnode->pop->event[msa_index], gnode->event);

        gnode->pop->event_count[msa_index]++;

        snode->seqin_count[msa_index]++;
      }
      
    }
  }
  return changed_count > 0;
}

static double lnprior_species_model(stree_t * stree)
{
  double p;

  switch (opt_delimit_prior)
  {
    case BPP_DELIMIT_PRIOR_DIRICHLET:
      p = 1.0;
      break;

    case BPP_DELIMIT_PRIOR_UNIFORM:
      p = 1.0 / histories(stree);
      break;

    default:
      fatal("Unknown species delimitation prior");
  }

  p = (p < 1e-300) ? -500 : log(p);
  
  return p;
}

long prop_split(gtree_t ** gtree,
                stree_t * stree,
                locus_t ** locus,
                double pr_split,
                long * param_count)
{
  int fsplit_count = 0;
  int fjoin_count = 0;
  int tau_count = 0;
  unsigned int i,k;
  long accepted = 0;
  double tau_new;
  double tau_upper;
  double pbetatau = 2;
  double qbetatau = 8;
  double thetafactor = 1;

  /* 1. Initialize lnacceptance */
  double lnacceptance = log((1-pr_split)/pr_split);

#if 0
  int * feasible = (int *)xmalloc((stree->tip_count + stree->inner_count) *
                           sizeof(int));
#endif

  double oldprior = lnprior_species_model(stree);

  /* 2. Count feasible nodes for splitting at source species tree */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau > 0)
      tau_count++;
    else
    {
      if (stree->nodes[i] == stree->root  || stree->nodes[i]->parent->tau > 0)
        feasible[fsplit_count++] = i;
    }
  }

  if (!fsplit_count)
    return 2;

//  printf("roottau = %f\n", stree->root_age);
    
  /* 3. Randomly select a feasible node to split, and compute an upper age */
  double r = legacy_rndu();
//  printf("random: %f\n", r);
  snode_t * node = stree->nodes[feasible[(int)(fsplit_count*r)]];
  if (node == stree->root)
    tau_upper = stree->root_age * 0.6;
  else
    tau_upper = node->parent->tau;

  /* 4. Change the age of the node, and update lnacceptance */
  node->old_tau = node->tau;
  node->tau = tau_new = tau_upper * legacy_rndbeta(pbetatau,qbetatau);
  lnacceptance -= log_pdfbeta(tau_new,pbetatau,qbetatau,tau_upper);

  /* TODO: No theta */

  /* 5. Now update population sizes for the two children, and then update
     lnacceptance */

  /* Store the left child theta, and update it according to RJ algorithm */
  node->left->old_theta = node->left->theta;

  if (!opt_rjalgorithm)
  {
    node->left->theta = node->theta*exp(opt_finetune_rja*(legacy_rndu() - 0.5));
    thetafactor *= opt_finetune_rja * node->left->theta;
  }
  else
  {
    node->left->theta = legacy_rndgamma(opt_finetune_rja) /
                        (opt_finetune_rja/(opt_finetune_rjm*node->theta));
    thetafactor /= pdf_gamma(node->left->theta,
                             opt_finetune_rja,
                             opt_finetune_rja/(opt_finetune_rjm*node->theta));
  }

  lnacceptance += log_pdfinvgamma(node->left->theta,
                                  opt_theta_alpha,
                                  opt_theta_beta);

  /* Store the right child theta, and update it according to RJ algorithm */
  node->right->old_theta = node->right->theta;

  if (!opt_rjalgorithm)
  {
    node->right->theta = node->theta*exp(opt_finetune_rja*(legacy_rndu() - 0.5));
    thetafactor *= opt_finetune_rja * node->right->theta;
  }
  else
  {
    node->right->theta = legacy_rndgamma(opt_finetune_rja) /
                         (opt_finetune_rja/(opt_finetune_rjm*node->theta));
    thetafactor /= pdf_gamma(node->right->theta,
                             opt_finetune_rja,
                             opt_finetune_rja/(opt_finetune_rjm*node->theta));
  }

  lnacceptance += log_pdfinvgamma(node->right->theta,
                                  opt_theta_alpha,
                                  opt_theta_beta);

  /* Update lnacceptance with new species tree prior */
  lnacceptance += lnprior_species_model(stree) - oldprior;

  if (node == stree->root)
    lnacceptance += log_pdfinvgamma(tau_new,opt_tau_alpha,opt_tau_beta);
  else
    lnacceptance += log(tau_count/stree->root->tau);    /* Eq 2 in YR2010 */

  /* 6. Count number of feasible nodes for joining at target */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * u = stree->nodes[i];
    if (u->tau== 0) continue;

    if ((!u->left->left || u->left->tau== 0) &&
        (!u->right->left || u->right->tau == 0))
      fjoin_count++;
  }

  lnacceptance += log((double)fsplit_count/fjoin_count * thetafactor);

  /* allocate storage for partials */

  /* 7. Call rubber-band to update all children */
  for (i = 0; i < stree->locus_count; ++i)
  {
    int j = rubber_proportional(stree,
                                node,
                                gtree,
                                tau_upper,
                                0,
                                tau_new,
                                i,
                                &lnacceptance);

    gtree[i]->old_logl = gtree[i]->logl;
    if (j)
    {
      locus_update_matrices_jc69(locus[i],
                                 nodevec+nodevec_offset[i],
                                 nodevec_count[i]);

      /* TODO: Never call functions like propose_age that change travbuffer
         from gtree.c */
      partials[i] = gtree_return_partials(gtree[i]->root,
                                          i,
                                          partials_count+i);
      for (k = 0; k < partials_count[i]; ++k)
        partials[i][k]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                   partials[i][k]->clv_index);

      /* update partials */
      locus_update_partials(locus[i],partials[i],partials_count[i]);

      
      /* evaluate log-likelihood */
      unsigned int param_indices[1] = {0};
      double logl = locus_root_loglikelihood(locus[i],
                                             gtree[i]->root,
                                             param_indices,
                                             NULL);
      gtree[i]->logl = logl;
    }
    else
    {
    //  assert(0);
    }

    double logpr = gtree[i]->logpr;

#if 1
    /* update log-pr */
    logpr -= node->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node,i);
    logpr -= node->left->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node->left,i);
    logpr -= node->right->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node->right,i);
#else
    logpr = gtree_logprob(stree,i);
#endif
    gtree[i]->old_logpr = gtree[i]->logpr;
    gtree[i]->logpr = logpr;

    lnacceptance += gtree[i]->logl  - gtree[i]->old_logl;
    lnacceptance += gtree[i]->logpr - gtree[i]->old_logpr;
  }

  if (opt_debug)
    printf("[Debug] (split) lnacceptance = %f\n", lnacceptance);

  if (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance))
  {
    /* accept */

    long dmodel_index = delimit_getindex(stree);

    delimit_setindex(dmodel_index);

    /* TODO: NoTheta option */
    *param_count += 1 + (node->left->theta > 0) + (node->right->theta > 0);


    accepted = 1;
    for (i = 0; i < stree->locus_count; ++i)
      for (k = 0; k < gtree[i]->tip_count + gtree[i]->inner_count; ++k)
      {
        gnode_t * tmp = gtree[i]->nodes[k];
        tmp->mark = 0;
      }
  }
  else
  {
    /* reject */

    node->tau = node->old_tau;
    node->left->theta  = node->left->old_theta;
    node->right->theta = node->right->old_theta;

    for (i = 0; i < stree->locus_count; ++i)
    {
      gtree[i]->logl  = gtree[i]->old_logl;
      gtree[i]->logpr = gtree[i]->old_logpr;

      /* TODO: Perhaps only nodevec needs to be traversed and not all nodes */
      for (k = 0; k < gtree[i]->tip_count + gtree[i]->inner_count; ++k)
      {
        gnode_t * tmp = gtree[i]->nodes[k];

        if (tmp->mark & MARK_AGE_UPDATE)
          tmp->time = tmp->old_time;

        if (tmp->mark & MARK_POP_CHANGE)
        {
          unlink_event(tmp,i);
          tmp->pop->event_count[i]--;

          tmp->pop = node;

          dlist_item_append(tmp->pop->event[i], tmp->event); /* equiv to snode->event[i] */

          tmp->pop->event_count[i]++;

          tmp->pop->seqin_count[i]++;
        }

        tmp->mark = 0;
      }
          
      locus_update_matrices_jc69(locus[i],
                                 nodevec+nodevec_offset[i],
                                 nodevec_count[i]);

      for (k = 0; k < partials_count[i]; ++k)
        partials[i][k]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                   partials[i][k]->clv_index);

      node->logpr_contrib[i] = node->old_logpr_contrib[i];
      node->left->logpr_contrib[i] = node->left->old_logpr_contrib[i];
      node->right->logpr_contrib[i] = node->right->old_logpr_contrib[i];
    }
  }

  #if 0
  /* CROSS VALIDATION */

  for (i = 0; i < stree->locus_count; ++i)
  {
    unsigned int j = 0;
    gtree_t * gt = gtree[i];

    gnode_t ** gt_nodes = (gnode_t **)xmalloc((gt->tip_count + gt->inner_count)*
                                              sizeof(gnode_t *));
    k=0;
    for (j = 0; j < gt->tip_count + gt->inner_count; ++j)
      if (gt->nodes[j]->parent)
        gt_nodes[k++] = gt->nodes[j];
    locus_update_matrices_jc69(locus[i],gt_nodes,k);

    gtree_all_partials(gt->root,gt_nodes,&k);
    for (j = 0; j < k; ++j)
      gt_nodes[j]->clv_index = SWAP_CLV_INDEX(gt->tip_count,gt_nodes[j]->clv_index);

    locus_update_partials(locus[i],gt_nodes,k);

    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    double logl = locus_root_loglikelihood(locus[i],gt->root,param_indices,NULL);


    double logpr = gtree_logprob(stree,i);
    printf("VERIFICATION %d logl: %f logpr: %f\n", i, logl, logpr);
    free(gt_nodes);
  }
  #endif

  for (i = 0; i < stree->locus_count; ++i)
    nodevec_count[i] = 0;

  return accepted;
}

long prop_join(gtree_t ** gtree,
               stree_t * stree,
               locus_t ** locus,
               double pr_split,
               long * param_count)
{
  int fsplit_count = 0;
  int fjoin_count = 0;
  int tau_count = 0;
  long accepted = 0;
  double tau_upper;
  double pbetatau = 2;
  double qbetatau = 8;
  unsigned int i,k;
  double thetafactor = 1;

  dbg_rj_join++;
  /* 1. Initialize lnacceptance */
  double lnacceptance = log((1-pr_split)/pr_split);

  double oldprior = lnprior_species_model(stree);

  /* 2. Count feasible nodes for joining at source species tree. A node is
        feasible for joining if each of its children is either a tip or already
        joined  */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * u = stree->nodes[i];

    if (u->tau == 0) continue;

    tau_count++;

    if ((!u->left->left || u->left->tau== 0) &&
        (!u->right->left || u->right->tau == 0))
      feasible[fjoin_count++] = i;
  }

  if (!fjoin_count)
    return 2;

  /* 3. Randomly select a feasible node to join, and compute an upper age */
  double r = legacy_rndu();
  snode_t * node = stree->nodes[feasible[(int)(fjoin_count*r)]];
  if (node == stree->root)
    tau_upper = stree->root_age * 0.6;
  else
    tau_upper = node->parent->tau;

  if (node->tau >= tau_upper) return 2;

  /* 4. Change the age of the node, and update lnacceptance */
  lnacceptance += log_pdfbeta(node->tau,pbetatau,qbetatau,tau_upper);

  /* TODO: No theta */

  /* Store the left child theta, and update it according to RJ algorithm */
  node->left->old_theta = node->left->theta;

  if (node->left->theta > 0)
  {
    if (!opt_rjalgorithm)
    {
      double y = exp(opt_finetune_rja*0.5);
      if (node->left->theta < node->theta/y ||
          node->left->theta > node->theta*y)
        return 2;  /* move disallowed */
      thetafactor /= opt_finetune_rja * node->left->theta;
    }
    else
    {
      thetafactor *= pdf_gamma(node->left->theta,
                               opt_finetune_rja,
                               opt_finetune_rja/(opt_finetune_rjm*node->theta));
    }
    lnacceptance -= log_pdfinvgamma(node->left->theta,
                                    opt_theta_alpha,
                                    opt_theta_beta);
  }

  /* Now store the right child theta, and update it according to RJ algorithm */
  node->right->old_theta = node->right->theta;

  if (node->right->theta > 0)
  {
    if (!opt_rjalgorithm)
    {
      double y = exp(opt_finetune_rja*0.5);
      if (node->right->theta < node->theta/y ||
          node->right->theta > node->theta*y)
        return 2;  /* move disallowed */
      thetafactor /= opt_finetune_rja * node->right->theta;
    }
    else
    {
      thetafactor *= pdf_gamma(node->right->theta,
                               opt_finetune_rja,
                               opt_finetune_rja/(opt_finetune_rjm*node->theta));
    }
    lnacceptance -= log_pdfinvgamma(node->right->theta,
                                    opt_theta_alpha,
                                    opt_theta_beta);
  }

  /* save old taus in case of rejection */
  node->old_tau = node->tau;
  node->left->old_tau = node->left->tau;
  node->right->old_tau = node->right->tau;

  node->tau = 0;
  node->left->theta = node->right->theta = -1;

  /* Update lnacceptance with new species tree prior */
  lnacceptance += lnprior_species_model(stree) - oldprior;

  if (node == stree->root)
    lnacceptance -= log_pdfinvgamma(node->old_tau,opt_tau_alpha,opt_tau_beta);
  else
    lnacceptance -= log((tau_count-1)/stree->root->tau);    /* Eq 2 in YR2010 */
  
  /* count feasible nodes for splitting at target. A node is feasible if its
     father is already split. */
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * u = stree->nodes[i];

    if (u->tau == 0 && (u == stree->root || u->parent->tau> 0))
      fsplit_count++;
  }

  lnacceptance += log((double)fjoin_count/fsplit_count * thetafactor);

  /* allocate storage for partials */

  /* 7. Call rubber-band to update all children */
  for (i = 0; i < stree->locus_count; ++i)
  {
    int j = rubber_proportional(stree,
                                node,
                                gtree,
                                tau_upper,
                                node->old_tau,
                                0,
                                i,
                                &lnacceptance);

    gtree[i]->old_logl = gtree[i]->logl;
    if (j)
    {
      locus_update_matrices_jc69(locus[i],
                                 nodevec+nodevec_offset[i],
                                 nodevec_count[i]);

      /* TODO: Never call functions like propose_age that change travbuffer
         from gtree.c */
      partials[i] = gtree_return_partials(gtree[i]->root,
                                          i,
                                          partials_count+i);

      for (k = 0; k < partials_count[i]; ++k)
        partials[i][k]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                   partials[i][k]->clv_index);

      /* update partials */
      locus_update_partials(locus[i],partials[i],partials_count[i]);

      
      /* evaluate log-likelihood */
      unsigned int param_indices[1] = {0};
      double logl = locus_root_loglikelihood(locus[i],
                                             gtree[i]->root,
                                             param_indices,
                                             NULL);
      gtree[i]->logl = logl;
    }
    else
    {
      
      //assert(0);
    }

    double logpr = gtree[i]->logpr;

#if 1
    /* update log-pr */
    logpr -= node->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node,i);
    logpr -= node->left->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node->left,i);
    logpr -= node->right->logpr_contrib[i];
    logpr += gtree_update_logprob_contrib(node->right,i);
#else
    logpr = gtree_logprob(stree,i);
#endif
    gtree[i]->old_logpr = gtree[i]->logpr;
    gtree[i]->logpr = logpr;

    lnacceptance += gtree[i]->logl  - gtree[i]->old_logl;
    lnacceptance += gtree[i]->logpr - gtree[i]->old_logpr;
  }

  if (opt_debug)
    printf("[Debug] (join) lnacceptance = %f\n", lnacceptance);

  if (lnacceptance >= 0 || legacy_rndu() < exp(lnacceptance))
  {
    /* accepted */

    long dmodel_index = delimit_getindex(stree);

    delimit_setindex(dmodel_index);

    accepted = 1;

    /* TODO: No theta option */
    *param_count -= 1 + (node->left->old_theta>0) + (node->right->old_theta>0);

    /* reset mark */
    for (i = 0; i < stree->locus_count; ++i)
      for (k = 0; k < gtree[i]->tip_count + gtree[i]->inner_count; ++k)
        gtree[i]->nodes[k]->mark = 0;

  }
  else
  {
    /* rejected */

    /* restore old tau */
    node->tau = node->old_tau;
    node->left->tau = node->left->old_tau;
    node->right->tau = node->right->old_tau;
    node->left->theta = node->left->old_theta;
    node->right->theta = node->right->old_theta;

    for (i = 0; i < stree->locus_count; ++i)
    {
      for (k = 0; k < gtree[i]->tip_count + gtree[i]->inner_count; ++k)
      {
        gnode_t * tmp = gtree[i]->nodes[k];

        if (tmp->mark & MARK_AGE_UPDATE)
          tmp->time = tmp->old_time; 

        if (tmp->mark & MARK_POP_CHANGE)
        {
          unlink_event(tmp,i);
          tmp->pop->event_count[i]--;
          tmp->pop->seqin_count[i]--;

          tmp->pop = tmp->old_pop;

          dlist_item_append(tmp->pop->event[i], tmp->event); /* equiv to snode->event[i] */

          tmp->pop->event_count[i]++;
        }

        /* reset marks */
        tmp->mark = 0;
      }

      /* restore logl and logpr for each gene tree */
      if (gtree[i]->logl != gtree[i]->old_logl)
      {
        locus_update_matrices_jc69(locus[i],
                                   nodevec+nodevec_offset[i],
                                   nodevec_count[i]);

        for (k = 0; k < partials_count[i]; ++k)
          partials[i][k]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                     partials[i][k]->clv_index);

        gtree[i]->logl = gtree[i]->old_logl;
      }
      gtree[i]->logpr = gtree[i]->old_logpr;

      node->logpr_contrib[i] = node->old_logpr_contrib[i];
      node->left->logpr_contrib[i] = node->left->old_logpr_contrib[i];
      node->right->logpr_contrib[i] = node->right->old_logpr_contrib[i];
    }
  }
  for (i = 0; i < stree->locus_count; ++i)
    nodevec_count[i] = 0;

  return accepted;
}
