### test_preprocessing_methods.pl ##################################################################

### NOTES ##########################################################################################
# Copied from ~/svn/Collaborators/RobBristow/nanostring_validation/normalization/

### PREAMBLE #######################################################################################
use warnings;
use strict;
#use BoutrosLab::Utilities::SGE::JobGroup;
use HPCI;
use File::Temp;
use File::pushd;
use Test::More ; # tests => ##
use Test::Exception;
use Test::File::Contents;
use MooseX::Types::Path::Class qw(Dir File);
use File::ShareDir;
use FindBin;
#use File::ShareDir

# get args
my $group_name = $ARGV[0];
my $matched = $ARGV[1];
my $cnas = $ARGV[2];

my $log_dir = "/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV/logs/";
my $script_dir = "/u/dsendorek/svn/Resources/code/R/NanoStringNormCNV/trunk/NanoStringNormCNV/inst/";
my $script_file = $script_dir . "assess_preprocessing_methods_2.0.R";

my $sge = HPCI->group(
	cluster => 'SGE',
	name => $group_name,
	base_dir => $log_dir,
	max_concurrent => 250
	);
#my $sge = BoutrosLab::Utilities::SGE::JobGroup->new(collect_job_stats => 1);

my @modules = ('R-BL');
my $num_jobs = 0;

my $vis = 0;

# my ($perchip, $ccn, $bc, $scc, $inv, $oth, $col) = (0,0,0,0,0,0,0);
for (my $perchip = 0; $perchip <= 1; ++$perchip) {
	for(my $ccn = 0; $ccn <= 2; ++$ccn){
		for(my $bc = 0; $bc <= 3; ++$bc){
			for(my $scc = 0; $scc <= 4; ++$scc){
				for(my $inv = 0; $inv <= 1; ++$inv){
					# for(my $matched = 0; $matched <= 1; ++$matched){
						for(my $oth = 0; $oth <= 3; ++$oth){
							# for(my $cnas = 0; $cnas <= 3; ++$cnas){
								for(my $col = 0; $col <= 1; ++$col){

									# skipping because this is the same as 'cnas' set to 1
									if ($matched == 1 && $cnas == 0) {
										next;
										}
									
									# stage job
									my $job_name = "perchip${perchip}_ccn${ccn}_bc${bc}_scc${scc}_inv${inv}_oth${oth}_matched${matched}_cnas${cnas}_col${col}_vis${vis}";
									print $job_name . "\n";
									$sge->stage(
										command => "Rscript " . $script_file . " --perchip $perchip --ccn $ccn --bc $bc --scc $scc --inv $inv --oth $oth --matched $matched --cnas $cnas --col $col --vis $vis",
										name => $job_name,
										modules_to_load => \@modules,
										should_save_script => 0,
										resources_required => {h_vmem => '4G'}
										);
									++$num_jobs;

									}
								}
						# 	}
						# }
					}
				}
			}
		}
	}

# launch everything
print "Loaded ${num_jobs} jobs\n";
my %res = %{$sge->execute()};

# check execution
print("Should have finished executing\n");
open(my $status_file, '>' . $log_dir . 'failed_jobs_' . $group_name . '.txt');
my $ret_val = 0;
foreach my $name (keys %res) {
	my $stage = $res{$name};
	my $stage_length = scalar @{$stage};
	if ($stage_length == 0 || $stage->[$stage_length - 1]->{'exit_status'} !~ m/^\s*0\s*$/g) {
		print $status_file "Failure of stage $name\n";
		$ret_val = 1;
		}
	}
