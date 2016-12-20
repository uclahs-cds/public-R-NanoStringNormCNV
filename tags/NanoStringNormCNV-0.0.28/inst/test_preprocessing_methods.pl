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

my $sge = HPCI->group(cluster => 'SGE', logdir => 'logs/', name => 'parameterize_ns');
#my $sge = BoutrosLab::Utilities::SGE::JobGroup->new(collect_job_stats => 1);

my @modules = ('R-BL');
my $num_jobs = 0;
for (my $perchip = 0; $perchip <= 1; ++$perchip) {
	for(my $ccn = 0; $ccn <= 2; ++$ccn){
		for(my $bc = 0; $bc <= 3; ++$bc){
			for(my $scc = 0; $scc <= 4; ++$scc){
				for(my $inv = 0; $inv <= 1; ++$inv){
					for(my $matched = 0; $matched <= 1; ++$matched){
						for(my $oth = 0; $oth <= 3; ++$oth){
							for(my $cnas = 0; $cnas <= 3; ++$cnas){
								for(my $col = 0; $col <= 1; ++$col){
									my $job_name = "perchip${perchip}_ccn${ccn}_bc${bc}_scc${scc}_inv${inv}_oth${oth}_matched${matched}_cnas${cnas}_col${col}";
									print $job_name . "\n";
									$sge->stage(
										command => "Rscript assess_preprocessing_methods.R --perchip $perchip --ccn $ccn --bc $bc --scc $scc --inv $inv --oth $oth --matched $matched --cnas $cnas --col $col",
										name => $job_name,
										modules_to_load => \@modules,
										should_save_script => 0,
										resources_required => {h_vmem => '4G'}
										);
									++$num_jobs;
									}
								}
							}
						}
					}
				}
			}
		}
	}

print "Loaded ${num_jobs} jobs\n";

my %res = %{$sge->execute()};
print("Should have finished executing");
open(my $status_file, '>failed_jobs.txt');
my $ret_val = 0;
foreach my $name (keys %res) {
	my $stage = $res{$name};
	my $stage_length = scalar @{$stage};
	if ($stage_length == 0 || $stage->[$stage_length - 1]->{'exit_status'} !~ m/^\s*0\s*$/g) {
		print $status_file "Failure of stage $name\n";
		$ret_val = 1;
		}
	}






