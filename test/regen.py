#!/usr/bin/env python

# Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
# Department of Genetics, Evolution and Environment,
# University College London, Gower Street, London WC1E 6BT, England

from subprocess import Popen, PIPE, call

import sys, stat, os
import tempfile
import time
import shutil
import urllib

#define hash for stable bpp version

opt_stable_hash = "6caaa426d0a2641722c6f996d5ccb219d1987542"

#define path to BPP binary
opt_bpp_bin = ""

# define each testbed and 

opt_testsuite_small_desc = "small"
opt_testsuite_small = [                 # [path-to-test,description]
   ["testbed/small/1",  "small-A00-1"],
   ["testbed/small/2",  "small-A00-2"],
   ["testbed/small/3",  "small-A00-3"],
   ["testbed/small/4",  "small-A00-4"],
   ["testbed/small/5",  "small-A00-5"],
   ["testbed/small/6",  "small-A00-6"],
   ["testbed/small/7",  "small-A00-7"],
   ["testbed/small/8",  "small-A00-8"],
   ["testbed/small/9",  "small-A00-9"],
   ["testbed/small/10", "small-A00-10"],
   ["testbed/small/11", "small-A00-11"],
   ["testbed/small/12", "small-A00-12"],
   ["testbed/small/13", "small-A00-13"],
   ["testbed/small/14", "small-A00-14"],
   ["testbed/small/15", "small-A00-15"],
   ["testbed/small/16", "small-A00-16"],
   ["testbed/small/17", "small-A01-17"],
   ["testbed/small/18", "small-A01-18"],
   ["testbed/small/19", "small-A01-19"],
   ["testbed/small/20", "small-A01-20"],
   ["testbed/small/21", "small-A01-21"],
   ["testbed/small/22", "small-A01-22"],
   ["testbed/small/23", "small-A01-23"],
   ["testbed/small/24", "small-A01-24"],
   ["testbed/small/25", "small-A01-25"],
   ["testbed/small/26", "small-A01-26"],
   ["testbed/small/27", "small-A01-27"],
   ["testbed/small/28", "small-A01-28"],
   ["testbed/small/29", "small-A01-29"],
   ["testbed/small/30", "small-A01-30"],
   ["testbed/small/31", "small-A01-31"],
   ["testbed/small/32", "small-A01-32"],
   ["testbed/small/33", "small-A01-33"],
   ["testbed/small/34", "small-A01-34"],
   ["testbed/small/35", "small-A01-35"],
   ["testbed/small/36", "small-A01-36"],
   ["testbed/small/37", "small-A01-37"],
   ["testbed/small/38", "small-A01-38"],
   ["testbed/small/39", "small-A01-39"],
   ["testbed/small/40", "small-A01-40"],
   ["testbed/small/41", "small-A01-41"],
   ["testbed/small/42", "small-A01-42"],
   ["testbed/small/43", "small-A01-43"],
   ["testbed/small/44", "small-A01-44"],
   ["testbed/small/45", "small-A01-45"],
   ["testbed/small/46", "small-A01-46"],
   ["testbed/small/47", "small-A01-47"],
   ["testbed/small/48", "small-A01-48"],
   ["testbed/small/49", "small-A10-49"],
   ["testbed/small/50", "small-A10-50"],
   ["testbed/small/51", "small-A10-51"],
   ["testbed/small/52", "small-A10-52"],
   ["testbed/small/53", "small-A10-53"],
   ["testbed/small/54", "small-A10-54"],
   ["testbed/small/55", "small-A10-55"],
   ["testbed/small/56", "small-A10-56"],
   ["testbed/small/57", "small-A10-57"],
   ["testbed/small/58", "small-A10-58"],
   ["testbed/small/59", "small-A10-59"],
   ["testbed/small/60", "small-A10-60"],
   ["testbed/small/61", "small-A10-61"],
   ["testbed/small/62", "small-A10-62"],
   ["testbed/small/63", "small-A10-63"],
   ["testbed/small/64", "small-A10-64"],
   ["testbed/small/65", "small-A10-65"],
   ["testbed/small/66", "small-A10-66"],
   ["testbed/small/67", "small-A10-67"],
   ["testbed/small/68", "small-A10-68"],
   ["testbed/small/69", "small-A10-69"],
   ["testbed/small/70", "small-A10-70"],
   ["testbed/small/71", "small-A10-71"],
   ["testbed/small/72", "small-A10-72"],
   ["testbed/small/73", "small-A10-73"],
   ["testbed/small/74", "small-A10-74"],
   ["testbed/small/75", "small-A10-75"],
   ["testbed/small/76", "small-A10-76"],
   ["testbed/small/77", "small-A10-77"],
   ["testbed/small/78", "small-A10-78"],
   ["testbed/small/79", "small-A10-79"],
   ["testbed/small/80", "small-A10-80"],
   ["testbed/small/81", "small-A10-81"],
   ["testbed/small/82", "small-A10-82"],
   ["testbed/small/83", "small-A10-83"],
   ["testbed/small/84", "small-A10-84"],
   ["testbed/small/85", "small-A10-85"],
   ["testbed/small/86", "small-A10-86"],
   ["testbed/small/87", "small-A10-87"],
   ["testbed/small/88", "small-A10-88"],
   ["testbed/small/89", "small-A10-89"],
   ["testbed/small/90", "small-A10-90"],
   ["testbed/small/91", "small-A10-91"],
   ["testbed/small/92", "small-A10-92"],
   ["testbed/small/93", "small-A10-93"],
   ["testbed/small/94", "small-A10-94"],
   ["testbed/small/95", "small-A10-95"],
   ["testbed/small/96", "small-A10-96"],
   ["testbed/small/97", "small-A10-97"],
   ["testbed/small/98", "small-A10-98"],
   ["testbed/small/99", "small-A10-99"],
   ["testbed/small/100","small-A10-100"],
   ["testbed/small/101","small-A10-101"],
   ["testbed/small/102","small-A10-102"],
   ["testbed/small/103","small-A10-103"],
   ["testbed/small/104","small-A10-104"],
   ["testbed/small/105","small-A10-105"],
   ["testbed/small/106","small-A10-106"],
   ["testbed/small/107","small-A10-107"],
   ["testbed/small/108","small-A10-108"],
   ["testbed/small/109","small-A10-109"],
   ["testbed/small/110","small-A10-110"],
   ["testbed/small/111","small-A10-111"],
   ["testbed/small/112","small-A10-112"],
   ["testbed/small/113","small-A11-113"],
   ["testbed/small/114","small-A11-114"],
   ["testbed/small/115","small-A11-115"],
   ["testbed/small/116","small-A11-116"],
   ["testbed/small/117","small-A11-117"],
   ["testbed/small/118","small-A11-118"],
   ["testbed/small/119","small-A11-119"],
   ["testbed/small/120","small-A11-120"],
   ["testbed/small/121","small-A11-121"],
   ["testbed/small/122","small-A11-122"],
   ["testbed/small/123","small-A11-123"],
   ["testbed/small/124","small-A11-124"],
   ["testbed/small/125","small-A11-125"],
   ["testbed/small/126","small-A11-126"],
   ["testbed/small/127","small-A11-127"],
   ["testbed/small/128","small-A11-128"],
   ["testbed/small/129","small-A11-129"],
   ["testbed/small/130","small-A11-130"],
   ["testbed/small/131","small-A11-131"],
   ["testbed/small/132","small-A11-132"],
   ["testbed/small/133","small-A11-133"],
   ["testbed/small/134","small-A11-134"],
   ["testbed/small/135","small-A11-135"],
   ["testbed/small/136","small-A11-136"],
   ["testbed/small/137","small-A11-137"],
   ["testbed/small/138","small-A11-138"],
   ["testbed/small/139","small-A11-139"],
   ["testbed/small/140","small-A11-140"],
   ["testbed/small/141","small-A11-141"],
   ["testbed/small/142","small-A11-142"],
   ["testbed/small/143","small-A11-143"],
   ["testbed/small/144","small-A11-144"],
   ["testbed/small/145","small-A11-145"],
   ["testbed/small/146","small-A11-146"],
   ["testbed/small/147","small-A11-147"],
   ["testbed/small/148","small-A11-148"],
   ["testbed/small/149","small-A11-149"],
   ["testbed/small/150","small-A11-150"],
   ["testbed/small/151","small-A11-151"],
   ["testbed/small/152","small-A11-152"],
   ["testbed/small/153","small-A11-153"],
   ["testbed/small/154","small-A11-154"],
   ["testbed/small/155","small-A11-155"],
   ["testbed/small/156","small-A11-156"],
   ["testbed/small/157","small-A11-157"],
   ["testbed/small/158","small-A11-158"],
   ["testbed/small/159","small-A11-159"],
   ["testbed/small/160","small-A11-160"],
   ["testbed/small/161","small-A11-161"],
   ["testbed/small/162","small-A11-162"],
   ["testbed/small/163","small-A11-163"],
   ["testbed/small/164","small-A11-164"],
   ["testbed/small/165","small-A11-165"],
   ["testbed/small/166","small-A11-166"],
   ["testbed/small/167","small-A11-167"],
   ["testbed/small/168","small-A11-168"],
   ["testbed/small/169","small-A11-169"],
   ["testbed/small/170","small-A11-170"],
   ["testbed/small/171","small-A11-171"],
   ["testbed/small/172","small-A11-172"],
   ["testbed/small/173","small-A11-173"],
   ["testbed/small/174","small-A11-174"],
   ["testbed/small/175","small-A11-175"],
   ["testbed/small/176","small-A11-176"]
]

opt_testsuite_ziheng_desc = "Edge cases reported by Ziheng Yang"
opt_testsuite_ziheng = [                 # [path-to-test,description]
   ["testbed/ziheng/1",  "ziheng-1"],
   ["testbed/ziheng/2",  "ziheng-2"],
   ["testbed/ziheng/3",  "ziheng-3"],
   ["testbed/ziheng/4",  "ziheng-4"]
]
# define test collections

opt_testbeds = [
   [opt_testsuite_small,opt_testsuite_small_desc],
   [opt_testsuite_ziheng,opt_testsuite_ziheng_desc]
]

## define architectures to test

opt_testarch = ["CPU","SSE","AVX","AVX2"]
opt_regenarch = ["AVX"]


##############################
# DO NOT MODIFY FROM HERE ON #
##############################

colors = {
   "default"  : "",
   "-"        : "\x1b[00m",
   "red"      : "\x1b[31;1m",
   "green"    : "\x1b[32;1m",
   "yellow"   : "\x1b[33;1m",
   "blue"     : "\x1b[34;1m",
   "magenta"  : "\x1b[35;1m",
   "cyan"     : "\x1b[36;1m",
   "bluebg"   : "\x1b[44;1m",
   "yellowbg" : "\x1b[43;2m",
   "cyanbg"   : "\x1b[46;2m",
   "magbg"    : "\x1b[45;2m"
 }

def get_stable_bpp(tmpdir):
  baseurl = 'https://github.com/bpp/bpp/archive/'
  url = baseurl + opt_stable_hash + '.zip'
  urllib.urlretrieve(url,tmpdir + '/' + opt_stable_hash + '.zip')

def header():
  sys.stdout.write(" _                   _  _   \n"
                   "| |                 | || |  \n"
                   "| |__  _ __  _ __   | || |_ \n"
                   "| '_ \| '_ \| '_ \  |__   _|\n"
                   "| |_) | |_) | |_) |    | |  \n"
                   "|_.__/| .__/| .__/     |_|  \n"
                   "      | |   | |             \n"
                   "      |_|   |_|             \n"
                   "\n"
                   "bpp 4 testing framework\n\n");
  
def has_colors(stream):
  if not hasattr(stream,"isatty"):
    return False
  if not stream.isatty():
    return False
  try:
    import curses
    curses.setupterm()
    return curses.tigetnum("colors") > 2
  except:
    return False

def ansiprint(color,text,breakline=0):
  if colors[color] and has_colors:
    sys.stdout.write(colors[color] + text + "\x1b[00m")
  else:
    sys.stdout.write(text)
  if breakline:
    sys.stdout.write("\n")

def testf(curtest,numtest,t,desc,arch):
  
  # create output directory
  outdir = t + "/out";
  if not os.path.exists(outdir):
    os.makedirs(outdir)

  refdir = t + "/ref";
  if not os.path.exists(refdir):
    os.makedirs(refdir)

  # remove old files if they exist in outdir
  if os.path.exists(outdir + "/mcmc.txt"):
    os.remove(outdir + "/mcmc.txt")
  if os.path.exists(outdir + "/out.txt"):
    os.remove(outdir + "/out.txt")
  if os.path.exists(outdir + "/FigTree.tre"):
    os.remove(outdir + "/FigTree.tre")

  # remove old files if they exist in refdir
  if os.path.exists(refdir + "/mcmc.txt"):
    os.remove(refdir + "/mcmc.txt")
  if os.path.exists(refdir + "/out.txt"):
    os.remove(refdir + "/out.txt")
  if os.path.exists(refdir + "/FigTree.tre"):
    os.remove(refdir + "/FigTree.tre")

  ctl = t + "/data/bpp.ctl"
  cmd = opt_bpp_bin + " --cfile " + ctl + " --arch "  + arch + " 2>tmperr >tmp"

  now = time.strftime("  %H:%M:%S")

  tstart = time.time()

  p1 = Popen(cmd, shell=True)
  os.waitpid(p1.pid,0)

  os.rename(outdir + "/mcmc.txt", refdir + "/mcmc.txt")
  os.rename(outdir + "/out.txt", refdir + "/out.txt")
  if os.path.exists(outdir + "/FigTree.tre"):
    os.rename(outdir + "/FigTree.tre", refdir + "/FigTree.tre")
    

  tend = time.time()
  runtime = tend - tstart

  ansiprint("-", "{:>3}/{:<3} ".format(curtest,numtest) + now)
  ansiprint("cyan", " {:<39} ".format(desc))

  runtime = "%.2f" % runtime
  ansiprint("cyan", "{:<14} ".format(runtime))
  print

   
def genrefs():
  total = 0;
  for tb in opt_testbeds:
    testlist = tb[0]
    desc = tb[1]
    total += len(testlist) 
  print(" %d tests found" % total)
  print(" %d arch sets" % len(opt_testarch))

  for arch in opt_regenarch:
    ansiprint("bluebg", "{:<80}"
               .format(arch.rjust(40+len(arch)/2)), True)
    ansiprint("yellowbg", "{:<7}   {:<8} {:<39} {:<14}       "
               .format(" ","Start", "Test", "Time [s]"),True)

    current = 0;
    for tb in opt_testbeds:
      testlist = tb[0]
      desc = tb[1]
      for t in testlist:
        test = t[0]
        testdesc = t[1];
        current = current+1
        testf(current,total,test,testdesc,arch)

def create_reference_data():
  
  # create temporary folder
  tmpdir = tempfile.mkdtemp()
  print "Created temporary folder " + tmpdir

  # get stable bpp
  print "Downloading BPP with hash " + opt_stable_hash + " ..."
  get_stable_bpp(tmpdir)

  # unzip
  zipfile = tmpdir + '/' + opt_stable_hash + '.zip'
  print "Unzipping " + zipfile
  p1 = Popen("unzip " + zipfile + " -d " + tmpdir + " > /dev/null 2>&1",shell=True)
  os.waitpid(p1.pid,0)

  # make
  print "Building executable..."
  p1 = Popen("make -C " + tmpdir + '/bpp-' + opt_stable_hash + '/src > /dev/null 2>&1',shell=True)
  os.waitpid(p1.pid,0)

  # generate ref solutions
  print "Generating reference solutions..."
  global opt_bpp_bin 
  opt_bpp_bin = tmpdir + '/bpp-' + opt_stable_hash + '/src/bpp'
  genrefs()

  # delete temp dir and all its contents
  shutil.rmtree(tmpdir)



if __name__ == "__main__":
  
  header()

  #if not os.path.isfile(os.path.expandvars(opt_bpp_bin)):
  #  print "BPP binary not found. Please update variable 'opt_bpp_bin' (line 29)"

  create_reference_data()
  #genrefs()
