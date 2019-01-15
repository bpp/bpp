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
  int work;

  /* thread parameters */
  locus_t ** locus;
  gtree_t ** gtree;
  stree_t * stree;
  long locus_first;
  long locus_count;

  /* return values */
  long proposals;
  long accepted;
} thread_info_t;

static thread_info_t * ti;
static pthread_attr_t attr;

static void * threads_worker(void * vp)
{
  long t = (long)vp;
  thread_info_t * tip = ti + t;

  pthread_mutex_lock(&tip->mutex);

  tip->proposals = 0;
  tip->accepted  = 0;

  /* loop until signalled to quit */
  while (tip->work >= 0)
  {
    /* wait for work available */
    if (tip->work == 0)
      pthread_cond_wait(&tip->cond, &tip->mutex);

    if (tip->work > 0)
    {
      /* TODO: call work */
      switch (tip->work)
      {
        case THREAD_WORK_GTAGE:
          gtree_propose_ages_parallel(tip->locus,
                                      tip->gtree,
                                      tip->stree,
                                      tip->locus_first,
                                      tip->locus_count,
                                      t,
                                      &tip->proposals,
                                      &tip->accepted);
          break;
        case THREAD_WORK_GTSPR:
          gtree_propose_spr_parallel(tip->locus,
                                     tip->gtree,
                                     tip->stree,
                                     tip->locus_first,
                                     tip->locus_count,
                                     t,
                                     &tip->proposals,
                                     &tip->accepted);
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

void threads_init()
{
  long t;

  assert(opt_threads <= opt_locus_count);

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
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;
    tip->locus = NULL;
    tip->gtree = NULL;
    tip->stree = NULL;

    /* allocate loci for thread t */
    tip->locus_first = loci_start;
    tip->locus_count = loci_per_thread + (loci_remaining > 0 ? 1 : 0);

    if (loci_remaining)
      --loci_remaining;
    loci_start += tip->locus_count;

    printf("Starting thread %ld (loci %ld-%ld)\n", t, tip->locus_first, tip->locus_first+tip->locus_count);
    
    pthread_mutex_init(&tip->mutex, NULL);
    pthread_cond_init(&tip->cond, NULL);
    if (pthread_create(&tip->thread, &attr, threads_worker, (void *)(long)t))
      fatal("Cannot create thread");
  }
}

void threads_wakeup(int work_type,
                    locus_t ** locus,
                    gtree_t ** gtree,
                    stree_t * stree,
                    double * rc)
{
  long t; 

  /* dynamic load distribution */
  /* Currently we do not do any dynamic load distribution. We assign a static
     workload at initialization, which then never changes */


  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    pthread_mutex_lock(&tip->mutex);

    tip->work = work_type;
    tip->locus = locus;
    tip->gtree = gtree;
    tip->stree = stree;

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

  if (work_type == 1 || work_type == 2)
  {
    long proposals = 0;
    long accepted = 0;
    for (t = 0; t < opt_threads; ++t)
    {
      thread_info_t * tip = ti+t;
      proposals += tip->proposals;
      accepted  += tip->accepted;
    }
    if (!accepted)
      *rc = (double)0;

    *rc = (double)accepted/proposals;
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
