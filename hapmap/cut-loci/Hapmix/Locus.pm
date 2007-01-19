package Hapmix::Locus;
use strict;

# Create a Locus instance from line as defined by HAPMIXMAP data file
# "SNPid" "NumAlleles"    "DistanceinMb"
sub new(@) {
    my $class = shift;
    my $self = {};
    # Read the data.
    $self->{SNP_ID} = shift;
    $self->{NUM_ALLELES} = shift;
    $self->{DISTANCE_IN_MB} = shift;
    my $position = shift;
    # Using a nasty thing: ("#" + 0) == 0 in Perl
    $self->{POSITION} = ($position + 0) + ($self->{DISTANCE_IN_MB} + 0) * 10e+6;
    $self->{NUMBER} = shift;
    if (0) {
        print "Created a Locus with:";
        print " " . $self->{SNP_ID};
        print " " . $self->{NUM_ALLELES};
        print " " . $self->{DISTANCE_IN_MB};
        print " " . $self->{POSITION};
        print " " . $self->{NUMBER};
        print "\n";
    }
    bless($self, $class);
    return $self;
}

sub snp_id {
    my $self = shift;
    return $self->{SNP_ID};
}

sub num_alleles {
    my $self = shift;
    return $self->{NUM_ALLELES};
}

sub distance_in_mb {
    my $self = shift;
    return $self->{DISTANCE_IN_MB};
}

sub position {
    my $self = shift;
    return $self->{POSITION};
}

sub get_line {
    my $self = shift;
    my $previous = shift;
    my $distance;
    if (not $previous) {
        $distance = "#";
    } else {
        $distance = $self->position() - $previous->position();
        $distance /= 10e+6;
        # Keep consistent format with the original locus files.
        $distance = sprintf("%1.8f", $distance);
    }
    return join("\t", (
            $self->snp_id(),
            $self->num_alleles(),
            $distance));
}

1;
