#!/usr/bin/perl

# =============================================================================
# Include
# =============================================================================
use strict;

# =============================================================================
# Main part
# =============================================================================

# reading arguments
if ($ARGV[0] eq "--help") {
  print STDOUT <DATA>;
  exit;
}

my $mem = 8;
my $cpu = 1;
my $cpu_defined = 0;
my $time = "3:0:0";
my $array_file;
while (scalar(@ARGV) > 0 and $ARGV[0] =~ m/^-/g) {
  my $carg = shift(@ARGV);
  if ($carg =~ m/^-R$/g) {
    $mem = shift(@ARGV);
  }
  elsif ($carg =~ m/^-c$/g) { 
    $cpu = shift(@ARGV);
    $cpu_defined = 1;
  }
  elsif ($carg =~ m/^-t$/g) { 
    $time = shift(@ARGV);
  }
  elsif ($carg =~ m/^-a$/g) { 
    $array_file = shift(@ARGV);
  }
}

if ($cpu_defined == 0) {
  $cpu = int($mem/8);
}

print STDERR "** Running job with $mem GB of memory, $cpu CPUs **\n";



# Job array
my $array_length = 0;
if (defined $array_file) {
    $array_length = `cat $array_file | wc -l`;
    chomp $array_length;
}
$array_length = $array_length;



# script
my $command = join(" ", @ARGV);
my $pwd_dir = `pwd`;
chomp $pwd_dir;

open(FILE,">qsub_$$.sh") or die("Cannot create script file\n");
print FILE "#!/bin/csh \n\n";
print FILE "#SBATCH --mem $mem"."GB\n";
print FILE "#SBATCH --time $time\n";
if ($cpu > 0) {
  print FILE "#SBATCH -c $cpu\n";
}
print FILE "#SBATCH -o qsub_$$.out\n";
print FILE "#SBATCH -e qsub_$$.err\n";
print FILE "#SBATCH -N 1\n\n";
print FILE "source ~/.my.cshrc;\n";
print FILE "cd $pwd_dir;\n";
print FILE "$command;\n";
if ($array_length > 0) {
    print FILE "echo task \$SLURM_ARRAY_TASK_ID out of $array_length;\n";
    print FILE "`cat $array_file | body.pl \$SLURM_ARRAY_TASK_ID \$SLURM_ARRAY_TASK_ID`;\n";
}
close(FILE);
system("chmod +x qsub_$$.sh");



# command line
my $pwd = `pwd`;
chomp $pwd;

my $qsub_command = "sbatch";
if (defined $array_file) {
    $qsub_command = $qsub_command." -a 0-$array_length";
}
$qsub_command = $qsub_command." \"./qsub_$$.sh\"";

# sending command to queue
print STDERR "Submitting command to queue: $qsub_command\n";
system("$qsub_command");
print STDERR "Done.\n";



# =============================================================================
# Subroutines
# =============================================================================

# ------------------------------------------------------------------------
# Help message
# ------------------------------------------------------------------------
__DATA__

sbatch.pl [options] <command>

Runs a command in the queue.

 -R <size>    Require <size> GB memory for job (default = 8)
 -t <time>    Time until job is killed (default = 3:0:0, short queue) 
 -c <cpu>     Number of cpus required per task (default = memory_size/8)
 -a <file>    Run a job-array with commands as specified in <file>
