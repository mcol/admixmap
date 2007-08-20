=head1 NAME

Hapmix::HapmixData -- class representing HAPMIXMAP data.

=head1 SYNOPSIS

    $hm = Hapmix::HapmixData->new($data_base_name);
    
    If the data files are: 
        - foo_genotypes.txt
        - foo_loci.txt

    Then you should use "foo" as constructor argument.

=head1 DESCRIPTION

This class is intended for manipulation of HAPMIXMAP data. It reads
a set of two files and can operate on named individuals and loci. It's
capable of saving new data file from a list of loci to include. It will
recalculate loci distances, so the newly written file will have correct
loci distance values.

=cut

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

# Save file in fastPHASE format
sub write_fastphase($) {
    my $self = shift;
    my $file_name = shift;
    my $no_loci_count = shift;
    open(FAST_PHASE, ">$file_name");
    print FAST_PHASE $self->{GENOTYPES}->get_no_individuals() . "\n";
    if (not $no_loci_count) {
        print FAST_PHASE $self->{LOCI}->get_no_loci() . "\n";
    }
    my @fastphase_genotypes = $self->{GENOTYPES}->get_fastphase_lines();
    foreach my $line (@fastphase_genotypes) {
        print FAST_PHASE $line . "\n";
    }
    close(FAST_PHASE);
}

# The following line is necessary.
1;
