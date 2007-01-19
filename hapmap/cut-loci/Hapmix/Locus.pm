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
    my $distance = $self->{DISTANCE_IN_MB};
    my $position = shift;
    # This produces a warning when -w flag is used.
    # However, no other efficient methos is available.
    if ($distance == 0 && $distance ne "0")  {
        # Not a number
        $self->{POSITION} = 0;
    } else {
        $self->{POSITION} = $position + $self->{DISTANCE_IN_MB} * 10e+6;
    }
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

sub number {
    my $self = shift;
    return $self->{NUMBER};
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

sub within_range($$) {
    my $self = shift;
    my $offset = shift;
    my $other_locus = shift;
    my $diff = $self->position() - $other_locus->position();
    # if (abs($diff) <= $offset) {
    #     print "me: " . $self->position() . ", diff $diff, offset: $offset\n";
    # }
    return (abs($diff) <= $offset);
}

1;
