/*
    Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

typedef struct thread_info_s
{
  pthread_t thread;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  
  /* type of work (0 no work) */
  volatile int work;

  /* thread parameters */
  long locus_first;
  long locus_count;

  thread_data_t td;

} thread_info_t;

static thread_info_t * ti;
static pthread_attr_t attr;

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
static void pin_to_core(long t)
{
  cpu_set_t cpuset;

  CPU_ZERO(&cpuset);
  CPU_SET(t,&cpuset);

  if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset))
    fatal("Error while pinning thread to core. "
          "Probably used more threads than available cores?");
}
#endif

static void * threads_worker(void * vp)
{
  long t = (long)vp;
  thread_info_t * tip = ti + t;

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
  pin_to_core((opt_threads_start-1)+(t*opt_threads_step));
#endif

  pthread_mutex_lock(&tip->mutex);

  /* loop until signalled to quit */
  while (tip->work >= 0)
  {
    /* wait for work available */
    if (tip->work == 0)
      pthread_cond_wait(&tip->cond, &tip->mutex);

    if (tip->work > 0)
    {
      /* work work! */
      switch (tip->work)
      {
        case THREAD_WORK_GTAGE:
          gtree_propose_ages_parallel(tip->td.locus,
                                      tip->td.gtree,
                                      tip->td.stree,
                                      tip->locus_first,
                                      tip->locus_count,
                                      t,
                                      &tip->td.proposals,
                                      &tip->td.accepted);
          break;
        case THREAD_WORK_GTSPR:
          gtree_propose_spr_parallel(tip->td.locus,
                                     tip->td.gtree,
                                     tip->td.stree,
                                     tip->locus_first,
                                     tip->locus_count,
                                     t,
                                     &tip->td.proposals,
                                     &tip->td.accepted);
          break;
        case THREAD_WORK_TAU:
          propose_tau_update_gtrees(tip->td.locus,
                                    tip->td.gtree,
                                    tip->td.stree,
                                    tip->td.snode,
                                    tip->td.oldage,
                                    tip->td.minage,
                                    tip->td.maxage,
                                    tip->td.minfactor,
                                    tip->td.maxfactor,
                                    tip->locus_first,
                                    tip->locus_count,
                                    tip->td.affected,
                                    tip->td.paffected_count,
                                    &tip->td.count_above,
                                    &tip->td.count_below,
                                    &tip->td.logl_diff,
                                    &tip->td.logpr_diff,
                                    t);
          break;
        case THREAD_WORK_MIXING:
          prop_mixing_update_gtrees(tip->td.locus,
                                    tip->td.gtree,
                                    tip->td.stree,
                                    tip->locus_first,
                                    tip->locus_count,
                                    tip->td.c,
                                    t,
                                    &tip->td.lnacceptance,
                                    NULL);
          break;
        case THREAD_WORK_ALPHA:
          locus_propose_alpha_parallel(tip->td.stree,
                                       tip->td.locus,
                                       tip->td.gtree,
                                       tip->locus_first,
                                       tip->locus_count,
                                       t,
                                       &tip->td.proposals,
                                       &tip->td.accepted);
          break;
        case THREAD_WORK_RATES:
          locus_propose_qrates_parallel(tip->td.stree,
                                        tip->td.locus,
                                        tip->td.gtree,
                                        tip->locus_first,
                                        tip->locus_count,
                                        t,
                                        &tip->td.proposals,
                                        &tip->td.accepted);
          break;
        case THREAD_WORK_FREQS:
          locus_propose_freqs_parallel(tip->td.stree,
                                       tip->td.locus,
                                       tip->td.gtree,
                                       tip->locus_first,
                                       tip->locus_count,
                                       t,
                                       &tip->td.proposals,
                                       &tip->td.accepted);
          break;
        case THREAD_WORK_BRATE:
          prop_branch_rates_parallel(tip->td.gtree,
                                     tip->td.stree,
                                     tip->td.locus,
                                     tip->locus_first,
                                     tip->locus_count,
                                     t,
                                     &tip->td.proposals,
                                     &tip->td.accepted);
          break;
        default:
          fatal("Unknown work function assigned to thread worker %ld", t);
                             
      }
        
      tip->work = 0;
      pthread_cond_signal(&tip->cond);
    }
  }
  pthread_mutex_unlock(&tip->mutex);

  pthread_exit(NULL);
}

void threads_pin_master()
{
  #if (defined(__linux__) && !defined(DISABLE_COREPIN))
    pin_to_core(opt_threads_start-1);
  #endif
}

void threads_init(locus_t ** locus)
{
  long i;
  long t;
  long patterns;

  assert(opt_threads <= opt_locus_count);

#if (defined(__linux__) && !defined(DISABLE_COREPIN))
  pin_to_core(opt_threads_start-1);
#endif

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* allocate memory for thread info */
  ti = (thread_info_t *)xmalloc((size_t)opt_threads * sizeof(thread_info_t));

  /* static load allocation */
  /* TODO: We discussed with Ziheng better ways to allocate loci on threads. One
     idea was to use some weights on the site count of each locus, i.e.
     sites^{2/3} as we need to account for the MSC density computation time and
     not only on the phylogenetic likelihood.

     For now we only distrubute an equal number of loci to each thread and
     cyclically distribute the overflow */

  long loci_per_thread = opt_locus_count / opt_threads;
  long loci_remaining = opt_locus_count % opt_threads;
  long loci_start = 0;


  /* init and create worker threads */
  printf("\nDistributing workload to threads:\n");
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;
    tip->td.locus = NULL;
    tip->td.gtree = NULL;
    tip->td.stree = NULL;

    /* allocate loci for thread t */
    tip->locus_first = loci_start;
    tip->locus_count = loci_per_thread + (loci_remaining > 0 ? 1 : 0);

    if (loci_remaining)
      --loci_remaining;
    loci_start += tip->locus_count;

    /* calculate number of site patterns send to current thread */
    patterns = 0;
    for (i = 0; i < tip->locus_count; ++i)
      patterns += locus[tip->locus_first+i]->sites;
      
    printf(" Thread %ld : loci [%ld-%ld), %ld patterns\n",
           t, tip->locus_first, tip->locus_first+tip->locus_count, patterns);
    
    pthread_mutex_init(&tip->mutex, NULL);
    pthread_cond_init(&tip->cond, NULL);
    if (pthread_create(&tip->thread, &attr, threads_worker, (void *)(long)t))
      fatal("Cannot create thread");
  }
}

void threads_wakeup(int work_type, thread_data_t * data)
{
  long t; 

  /* dynamic load distribution */
  /* Currently we do not do any dynamic load distribution. We assign a static
     workload at initialization, which then never changes */


  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    pthread_mutex_lock(&tip->mutex);

    memcpy(&tip->td,data,sizeof(thread_data_t));

    tip->work = work_type;

    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);
  }

  /* wait for threads to finish their work */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti+t;

    pthread_mutex_lock(&tip->mutex);
    while (tip->work > 0)
      pthread_cond_wait(&tip->cond,&tip->mutex);

    pthread_mutex_unlock(&tip->mutex);
  }

  if (work_type == THREAD_WORK_GTAGE ||
      work_type == THREAD_WORK_GTSPR ||
      work_type == THREAD_WORK_ALPHA ||
      work_type == THREAD_WORK_RATES ||
      work_type == THREAD_WORK_FREQS ||
      work_type == THREAD_WORK_BRATE)
  {
    long proposals = 0;
    long accepted = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;
      proposals += tip->td.proposals;
      accepted  += tip->td.accepted;
    }

    data->proposals = proposals;
    data->accepted = accepted;
  }
  else if (work_type == THREAD_WORK_TAU)
  {
    data->count_above = 0;
    data->count_below = 0;
    data->logl_diff = 0;
    data->logpr_diff = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;

      data->count_above += tip->td.count_above;
      data->count_below += tip->td.count_below;
      data->logl_diff   += tip->td.logl_diff;
      data->logpr_diff  += tip->td.logpr_diff;
    }
  }
  else if (work_type == THREAD_WORK_MIXING)
  {
    data->lnacceptance = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;
      data->lnacceptance += tip->td.lnacceptance;
    }
  }
  else
    assert(0);
}

void threads_exit()
{
  long t;

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    /* tell worker to quit */
    pthread_mutex_lock(&tip->mutex);
    tip->work = -1;
    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);

    /* wait for worker to quit */
    if (pthread_join(tip->thread, 0))
      fatal("Cannot join thread");

    pthread_cond_destroy(&tip->cond);
    pthread_mutex_destroy(&tip->mutex);
  }

  free(ti);
  pthread_attr_destroy(&attr);
}
