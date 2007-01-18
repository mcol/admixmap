package Hapmix::HapmixData;
use strict;
use Hapmix::Loci;
use Hapmix::Genotypes;
sub new {
    my $class = shift;
    my $self = {};
    # Read the data.
    my $data_base_name = shift;
    my $genotypes_file = "${data_base_name}_genotypes.txt";
    my $locus_file = "${data_base_name}_loci.txt";

    $self->{LOCI} = Hapmix::Loci->new($locus_file);
    $self->{GENOTYPES} = Hapmix::Genotypes->new($genotypes_file);

    bless($self, $class);
    return $self;
}

sub write_genotypes_by_loci_ids($\@) {
    my $self = shift;
    my $file_name = shift;
    my @ids = @_;
    $self->{GENOTYPES}->write_file_from_loci_ids($file_name, @ids);
}

# The following line is necessary.
1;
