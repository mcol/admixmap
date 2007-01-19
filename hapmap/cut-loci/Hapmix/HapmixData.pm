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
    # Ensure that they are in correct order.
    my @loci = map {$self->{LOCI}->get_locus_by_id($_) } @ids;
    @loci = sort { $a->position() <=> $b->position() } @loci;
    @ids = map { $_->snp_id() } @loci;
    $self->{GENOTYPES}->write_file_from_loci_ids($file_name, @ids);
}

sub write_loci_by_ids($\@) {
    my $self = shift;
    my $file_name = shift;
    my @ids = @_;
    $self->{LOCI}->write_file_from_loci_ids($file_name, @ids);
}

# Writes genotypes and loci file
sub write_by_loci_ids($\@) {
    my $self = shift;
    my $base_name = shift;
    my @names = @_;
    my $genotypes_file = "${base_name}_genotypes.txt";
    my $locus_file = "${base_name}_loci.txt";
    $self->write_genotypes_by_loci_ids($genotypes_file, @names);
    $self->write_loci_by_ids($locus_file, @names);
}

sub offset($@) {
    my $self = shift;
    my $offset = shift;
    my @ids = @_;
    return $self->{LOCI}->offset($offset, @ids);
}

# Consider given list a range definition
# Return a list of ids.
sub range_by_ids(@) {
    my $self = shift;
    my @ids = @_;
    my @minmax = $self->{LOCI}->get_min_max_loci(@ids);
    return $self->{LOCI}->get_range(@minmax);
}

# The following line is necessary.
1;
