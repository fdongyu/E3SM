#!/usr/bin/env python
# script for running E3SM-MMF simulations using the 2020 INICTE allocation (CLI115)
# Branch for this campaign: https://github.com/E3SM-Project/E3SM/tree/whannah/incite-2020
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os
import subprocess as sp
import datetime
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

# directory info
case_dir = os.getenv('HOME')+'/E3SM/Cases/'
src_dir  = os.getenv('HOME')+'/E3SM/E3SM_SRC1/'

# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True

stop_opt,stop_n,resub = 'ndays',2,0

# this stuff should be the same for all runs
compset = 'F-MMFXX'
arch    = 'GNUGPU'
ne,npg  = 45,2
grid    = f'ne{ne}pg{npg}_r05_oECv3'
crm_dx  = 3200 
crm_dt  = 10
rad_nx  = 4

timestamp = '20200930'

case = '.'.join(['INCITE2020',arch,grid,compset,'NLEV_50','CRMNX_32','CRMNY_32','MOMFB',timestamp]); num_nodes=1000
case = '.'.join(['INCITE2020',arch,grid,compset,'NLEV_50','CRMNX_32','MOMFB',timestamp]); num_nodes=30
# case = '.'.join(['E3SM',arch,grid,compset,'CRMNX_32','CRMNY_32','CRMDX_3200','CRMDT_10','NLEV_50','RADNX_4','MOMFB',init_date,'00'])

# Impose wall limits for Summits
if num_nodes>=  1: walltime =  '2:00'
if num_nodes>= 46: walltime =  '6:00'
if num_nodes>= 92: walltime = '12:00'
if num_nodes>=922: walltime = '24:00'

### add these modifiers to enable debug mode or state variable checks
# case = case+'.debug-on'
# case = case+'.checks-on'

# init_file_dir = '/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch/init_files'
init_file_dir = '/gpfs/alpine/scratch/hannah6/cli115/HICCUP/data/'
init_file_sst = f'HICCUP.sst_noaa.{init_date}.nc'
params = [p.split('_') for p in case.split('.')]
for p in params:
   if p[0]=='NLEV': init_file_atm = f'HICCUP.cami_mam3_Linoz_ne{ne}np4.L{p[1]}_alt.nc'

### land IC
# land_init_path = '/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch/init_files'
# land_init_file = 'CLM_spinup.ICRUCLM45.ne45pg2_r05_oECv3.20-yr.2010-10-01.clm2.r.2005-10-01-00000.nc'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

dtime = 20*60   # GCM physics time step

if 'dtime' in locals(): ncpl  = 86400 / dtime

if arch=='GNUCPU': task_per_node = 84
if arch=='GNUGPU': task_per_node = 12
#-------------------------------------------------------------------------------
# Define run command
#-------------------------------------------------------------------------------
# Set up terminal colors
class tcolor:
   ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,suppress_output=False,execute=True):
   if suppress_output : cmd = cmd + ' > /dev/null'
   msg = tcolor.GREEN + cmd + tcolor.ENDC
   print(f'\n{msg}')
   if execute: os.system(cmd)
   return
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   cmd = src_dir+'cime/scripts/create_newcase --case '+case_dir+case
   cmd = cmd + ' --compset '+compset+' --res '+grid
   cmd = cmd + ' --pecount '+str(num_nodes*tasks_per_node)+'x1 '
   if arch=='GNUCPU': cmd = cmd + ' -compiler gnu    -pecount '+str(num_nodes*task_per_node)+'x1 '
   if arch=='GNUGPU': cmd = cmd + ' -compiler gnugpu -pecount '+str(num_nodes*task_per_node)+'x1 '
   run_cmd(cmd)

   # Change run directory to be next to bld directory
   os.chdir(case_dir+case+'/')
   run_cmd('./xmlchange -file env_run.xml RUNDIR=\'$MEMBERWORK/cli115/e3sm_scratch/$CASE/run\' ' )
   
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE    -val {task_per_node} ')
   run_cmd(f'./xmlchange -file env_mach_pes.xml -id MAX_MPITASKS_PER_NODE -val {task_per_node} ')

#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   #-------------------------------------------------------
   if arch=='GNUGPU': 
      pcols = np.ceil( (ne**2*6*npg**2) / (num_nodes*task_per_node) )
      run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -pcols {int(pcols)} \" ' )
   #-------------------------------------------------------
   file = open('user_nl_eam','w')
   file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')
   file.close()
   #-------------------------------------------------------
   # Specify CRM and RAD columns
   params = [p.split('_') for p in case.split('.')]
   for p in params:
      if p[0]=='CRMNX': run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -crm_nx {p[1]} \" ')
      if p[0]=='CRMNY': run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -crm_ny {p[1]} \" ')
      if p[0]=='NLEV' and p[1] != '72' : 
         nlev = p[1]; crm_nz = None
         if nlev== '50': crm_nz =  '46'
         if nlev=='100': crm_nz =  '95'
         if nlev=='120': crm_nz = '115'
         if crm_nz is None: raise ValueError('No value of crm_nz specified')
         run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -nlev {p[1]} -crm_nz {crm_nz} \" ')
   # add in the values that are always specified above
   rad_ny = rad_nx if 'CRMNY' in case else 1
   run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -crm_dx {crm_dx} -crm_dt {cmr_dt}  \" ')
   run_cmd(f'./xmlchange --append -id EAM_CONFIG_OPTS -val \" -crm_nx_rad {rad_nx} -crm_ny_rad {rad_ny} \" ')
   #-------------------------------------------------------
   # Add special MMF options based on case name
   cpp_opt = ''
   if '.MOMFB.'  in case: cpp_opt += ' -DMMF_MOMENTUM_FEEDBACK'
   if 'debug-on' in case: cpp_opt += ' -DYAKL_DEBUG'

   if cpp_opt != '' :
      cmd  = f'./xmlchange --append -file env_build.xml -id EAM_CONFIG_OPTS'
      cmd += f' -val \" -cppdefs \'-DMMF_DIR_NS {cpp_opt} \'  \" '
      run_cmd(cmd)
   #-------------------------------------------------------
   # Set tasks and threads
   if arch == 'GNUGPU':
      atm_ntasks = task_per_node*num_nodes
      cmd = './xmlchange -file env_mach_pes.xml '
      cmd += f' NTASKS_ATM={atm_ntasks}'
      if ne==4:  alt_ntask = task_per_node*4
      if ne==45: alt_ntask = task_per_node*6*4
      cmd += f',NTASKS_LND={alt_ntask},NTASKS_CPL={alt_ntask}'
      cmd += f',NTASKS_OCN={alt_ntask},NTASKS_ICE={alt_ntask}'
      alt_ntask = task_per_node
      cmd += f',NTASKS_ROF={alt_ntask},NTASKS_WAV={alt_ntask},NTASKS_GLC={alt_ntask}'
      cmd += f',NTASKS_ESP=1,NTASKS_IAC=1'
      run_cmd(cmd)
      # run_cmd('./xmlchange -file env_mach_pes.xml NTHRDS_ATM=2,NTHRDS_CPL=2,NTHRDS_LND=1')
   #-------------------------------------------------------
   # 64_data format is needed for large numbers of columns (GCM or CRM)
   run_cmd('./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\" ')
   #-------------------------------------------------------
   # Run case setup
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')

#---------------------------------------------------------------------------------------------------
# Build
#---------------------------------------------------------------------------------------------------
if build : 
   # run_cmd('./xmlchange PIO_VERSION=2')

   if 'debug-on' in case : run_cmd('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 
   # Change inputdata from default due to permissions issue
   run_cmd('./xmlchange DIN_LOC_ROOT=/gpfs/alpine/cli115/scratch/hannah6/inputdata ')
   
   #-------------------------------------------------------
   # First query some stuff about the case
   #-------------------------------------------------------
   (din_loc_root   , err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
                                     stdout=sp.PIPE, shell=True).communicate()
   # (EAM_CONFIG_OPTS, err) = sp.Popen('./xmlquery EAM_CONFIG_OPTS -value', \
   #                                   stdout=sp.PIPE, shell=True).communicate()
   # EAM_CONFIG_OPTS = ' '.join(EAM_CONFIG_OPTS.split())   # remove extra spaces to simplify string query
   ntasks_atm = None
   (ntasks_atm     , err) = sp.Popen('./xmlquery NTASKS_ATM    -value', \
                                     stdout=sp.PIPE, shell=True, universal_newlines=True).communicate()
   ntasks_atm = float(ntasks_atm)
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   #------------------------------
   # Specify history output frequency and variables
   #------------------------------   
   file.write(' nhtfrq    = 0,-1,-6 \n') 
   file.write(' mfilt     = 1,24, 4 \n')     
   # if npg>0 : file.write(" fincl1    = 'DYN_T','DYN_Q','DYN_U' \n")
   file.write(" fincl2    = 'PS','PSL','TS'")
   file.write(             ",'PRECT','TMQ'")
   file.write(             ",'LHFLX','SHFLX'")             # surface fluxes
   file.write(             ",'FSNT','FLNT'")               # Net TOM heating rates
   file.write(             ",'FLNS','FSNS'")               # Surface rad for total column heating
   file.write(             ",'FSNTC','FLNTC'")             # clear sky heating rates for CRE
   file.write(             ",'LWCF','SWCF'")               # cloud radiative foricng
   file.write(             ",'TGCLDLWP','TGCLDIWP'")
   file.write(             ",'TAUX','TAUY'")                       # surface stress
   file.write(             ",'TBOT:I','QBOT:I','UBOT:I','VBOT:I'") # lowest model leve
   file.write(             ",'T900:I','Q900:I','U900:I','V900:I'") # 900mb data
   file.write(             ",'T850:I','Q850:I','U850:I','V850:I'") # 850mb data
   file.write(             ",'Z300:I','Z500:I'")
   file.write(             ",'OMEGA850:I','OMEGA500:I'")
   file.write('\n')
   file.write(" fincl3    = 'PS','T','Q','Z3'")            # 3D thermodynamic budget components
   file.write(             ",'U','V','OMEGA'")             # 3D velocity components
   file.write(             ",'CLOUD','CLDLIQ','CLDICE'")   # 3D cloud fields
   file.write(             ",'QRL','QRS'")                 # 3D radiative heating profiles
   file.write(             ",'MMF_TKE','MMF_TKEW','MMF_TKES' ")
   # if 'SP' in cld :
   #    file.write(         ",'SPDT','SPDQ','SPDQC','SPDQI'")               # CRM heating/moistening tendencies
   #    file.write(         ",'SPTLS','SPQTLS'")           # CRM large-scale forcing
   #    file.write(         ",'SPQPEVP'")                   # CRM rain evap 
   #    # file.write(         ",'SPMC'")                      # total mass flux
   #    # file.write(         ",'SPMCUP','SPMCDN'")           # CRM saturated mass fluxes
   #    # file.write(         ",'SPMCUUP','SPMCUDN'")         # CRM unsaturated mass fluxes
   #    file.write(         ",'SPTKE','SPTKES'")
   #    # file.write(         ",'CRM_WSX:I','CRM_WSY:I' ")
   
   #    file.write(         ",'CRM_QC:I','CRM_QI:I'")
   #    file.write(         ",'CRM_U:I','CRM_V:I','CRM_W:I' ")
   #    # if any(x in EAM_CONFIG_OPTS for x in ["SP_ESMT","SP_USE_ESMT","SPMOMTRANS"]) : 
   #    #    file.write(",'ZMMTU','ZMMTV','uten_Cu','vten_Cu' ")
   #    # if "SP_USE_ESMT" in EAM_CONFIG_OPTS : file.write(",'U_ESMT','V_ESMT'")
   #    # if "SPMOMTRANS"  in EAM_CONFIG_OPTS : file.write(",'UCONVMOM','VCONVMOM'")
   file.write('\n')
   #------------------------------
   # Other namelist stuff
   #------------------------------
   # file.write(' srf_flux_avg = 1 \n')              # sfc flux smoothing (for MMF stability)
   # file.write(f' crm_accel_factor = 4 \n')         # CRM acceleration factor (default is 2)

   if num_dyn<ntasks_atm:
      file.write(' dyn_npes = '+str(num_dyn)+' \n')   # limit dynamics tasks

   if 'checks-on' in case : file.write(' state_debug_checks = .true. \n')

   file.write(" inithist = \'ENDOFRUN\' \n") # ENDOFRUN / NONE

   if 'init_file_atm' in locals():
      file.write(f' ncdata = \'{init_file_dir}/{init_file_atm}\'\n')

   file.close()

   #-------------------------------------------------------
   # CLM namelist
   #-------------------------------------------------------
   if 'land_init_file' in locals():
      nfile = 'user_nl_elm'
      file = open(nfile,'w')
      file.write(f' finidat = \'{land_init_path}/{land_init_file}\' \n')
      file.close()
   
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   if 'ncpl' in locals(): run_cmd(f'./xmlchange ATM_NCPL={str(ncpl)}')
   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_QUEUE=batch,JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')

   # # Restart Frequency
   # run_cmd('./xmlchange -file env_run.xml REST_OPTION={stop_opt},REST_N={stop_n}')

   if continue_run :
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=TRUE ')   
   else:
      run_cmd('./xmlchange -file env_run.xml CONTINUE_RUN=FALSE ')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   run_cmd('./case.submit')

#---------------------------------------------------------------------------------------------------
# Print the case name again
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------