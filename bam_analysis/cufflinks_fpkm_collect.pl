#!/usr/bin/perl

# =============================================================================
# Include
# =============================================================================
use strict;

# =============================================================================
# Main part
# =============================================================================

# Reading arguments
if ($ARGV[0] eq "--help") {
  print STDOUT <DATA>;
  exit;
}

my $file_ref;
my $file_name = $ARGV[0];
if (length($file_name) < 1 or $file_name =~ /^-/) {
  $file_ref = \*STDIN;
}
else {
  shift(@ARGV);
  open(FILE, $file_name) or die("Could not open file '$file_name'.\n");
  $file_ref = \*FILE;
}

# parameters
my $t = 0;
my $org = "DAR";
my $n = 0;
my $x = 0;
while (scalar(@ARGV) > 0 and $ARGV[0] =~ m/^-/g) {
  my $carg = shift(@ARGV);
  if ($carg =~ m/^-t/g) {
    $t = 1;
  }
  if ($carg =~ m/^-o/g) {
    $org = shift(@ARGV);
  }
  elsif ($carg =~ m/^-n/g) {
    $n = 1;
  }
  elsif ($carg =~ m/^-x/g) {
    $x = 1;
  }
}

my $idx = "ENS".$org."X";
my $idt = "ENS".$org."T";
print STDERR ">> $idx, $idt\n";

my %fpkm_P;
my %fpkm_M;
my %fpkm_O;
while (<$file_ref>) {
  chomp $_;
  my ($tid,$gid,$gname,$fpkm) = split("\t",$_);
  my $id = "$gid\t$gname";

  if ($t > 0) {
    if ($tid =~ m/transcript:$idx/g) {
      $fpkm_P{$id}+=$fpkm;
    }
    elsif ($tid =~ m/transcript:$idt/g) {
      $fpkm_M{$id}+=$fpkm;
    }
    else {
      $fpkm_O{$id}+=$fpkm;
    }
  }
  elsif ($n > 0) {
    if ($tid =~ m/$gid:pre/g) {
      $fpkm_P{$id}+=$fpkm;
    }
    else {
      $fpkm_M{$id}+=$fpkm;
    }
  }
  elsif ($x > 0) {
    if ($tid =~ m/^pRNA\./g) {
      $fpkm_P{$id}+=$fpkm;
    }
    else {
      $fpkm_M{$id}+=$fpkm;
    }
  }
  else {
    if ($tid =~ m/transcript:pre/g) {
      $fpkm_P{$id}+=$fpkm;
    }
    elsif ($tid =~ m/transcript:mat/g) {
      $fpkm_M{$id}+=$fpkm;
    }
    else {
      $fpkm_O{$id}+=$fpkm;
    }
  }
}

my @id_list = keys(%fpkm_P);
push(@id_list, keys(%fpkm_M));
push(@id_list, keys(%fpkm_O));
my @uids = uniq(@id_list);

foreach my $i (@uids) {
  my $p = 0;
  my $m = 0;
  if (exists $fpkm_M{$i}) {
    $m = $fpkm_M{$i};
    if (exists $fpkm_P{$i}) {
      $p = $fpkm_P{$i};
    }
  }
  elsif (exists $fpkm_P{$i}) {
    $m = $fpkm_P{$i};
  }
  else {
    $m = $fpkm_O{$i};
  }
  print "$i\t$m\t$p\n";
}


# =============================================================================
# Subroutines
# =============================================================================
sub uniq(@) {
  my (@l) = @_;
  
  my %seen;
  foreach my $i (@l) {
    $seen{$i} = 1;
  }
  return keys(%seen);
}

# ------------------------------------------------------------------------
# Help message
# ------------------------------------------------------------------------
__DATA__

cufflinks_fpkm_collect.pl 

  Analyze isoform fpkms calculated by cufflinks, to produce
  pre-mRNA/mRNA estimates of expression.

  Input format: <transcript id> <gene id> <gene name> <FPKM>
  Output format: <gene id> <gene name> <mRNA:FPKM> <pre-mRNA:FPKM>
  
Options:
  -t       use ENS[organism]X/ENS[organism]T titles
  -o <id>  organism (default = DAR, zebrafish)
  -n       use :pre titles
  -x       use pRNA. titles
