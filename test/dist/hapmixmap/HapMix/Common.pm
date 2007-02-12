#!/usr/bin/perl
package HapMix::Common;

use strict;

sub new {
    my $class = shift;
    my $self = {};
    $self->{POPULATIONS} = {
        Eur => "CEU",
        Afr => "YRI",
        Asian => "JPT+HBT",
    };
    $self->{DIRECTORIES} = {
        source          => "01-source",
        pre_processed   => "02-pre-processed",
        training        => "03-training",
        testing         => "04-testing",
        post_processing => "05-post-processing",
    };
    $self->{SRC_FILE} = "source-genotypes.txt";
    $self->{HAPMAP_BASE_URL} = "http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/";
    bless($self, $class);
    return $self;
}

1;
