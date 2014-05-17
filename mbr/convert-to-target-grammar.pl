#!/usr/bin/perl -w
use strict;

# [X_16_25] ||| [X_16_22] um die [X_24_25] ||| [1] [2] ||| spectral=-1.371 DeletionRule=1.0 EgivenFCoherent=0.0 SampleCountF=0.0 CountEF=0.0 MaxLexFgivenE=0.0 MaxLexEgivenF=0.0 IsSingletonF=0.0 IsSingletonFE=0.0

while(<>) {
  chomp;
  my ($lhs, $sf, $se, $ss) = split / \|\|\| /;
  my @nts;
  my @ff = split /\s+/, $sf;
  for my $f (@ff) {
    if ($f =~ /\[[^]]+\]/) {
      push @nts, $f;
    }
  }
  my @ee = split /\s+/, $se;
  for (my $i=0; $i< scalar @ee; $i++) {
    if ($ee[$i] =~ /\[(\d+)\]/) {
      my $ind = $1 - 1;
      my $nt = $nts[$ind];
      $ee[$i] = $nt;
    }
  }
  my @scores = split /\s+/, $ss;
  my $spec = shift @scores;
  $spec =~ s/^spectral=//;
  $spec = exp($spec);
  print "$lhs ||| @ee ||| Marginal=$spec\n";
}


