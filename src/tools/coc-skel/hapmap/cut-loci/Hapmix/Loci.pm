=head1 NAME

Hapmix::Loci -- class representing a set of loci.

=head1 SYNOPSIS

    use Hapmap::Loci;

    $loci = Hapmap::Loci->new($locus_file);
    $loci->write_file_from_loci_ids($file_name, @loci_ids);

=head1 DESCRIPTION

This class is intended to be used from withing Hapmix::HapmixData
module.

=cut

package Hapmix::Loci;
use strict;
use Hapmix::Locus;
use List::Util qw(min max);

sub new {
    my $class = shift;
    my $self = {};
    # Read the data.
    my $locus_file = shift;

    open(LOCI, "<$locus_file")
        or die("Can't open the locus file: $locus_file.\n");
    my @lines = <LOCI>;
    close(LOCI);

    # Set the header line.
    # NB: Cannae do the following as $self not 'blessed' yet.
    # $self->header_line(shift(@lines));
    # Need to use $self->{SOMETHING}.
    $self->{HEADER_LINE} = shift(@lines);
    # my @headers = split(/\s+/, $self->{HEADER_LINE});

    # $self->{LOCI_ARRAY} = \@lines;
    $self->{LOCI} = {};
    $self->{POSITIONS} = {};
    $self->{NUMBERS} = {};
    my $position = 0;
    my $counter = 0;
    foreach my $line (@lines) {
        my @cols = split(" ", $line);
        my $locus = Hapmix::Locus->new(@cols, $position, $counter);
        $self->{LOCI}{$cols[0]} = $locus;
        $position = $locus->position();
        # Hash loci by positions
        $self->{POSITIONS}{$position} = $locus;
        $self->{NUMBERS}{$counter} = $locus;
        $counter++;
    }
    $self->{HIGHEST_NUMBER} = $counter - 1;
    undef @lines;
    bless($self, $class);
    # print "Loci constructed from $locus_file.\n";
    return $self;
}

sub get_lines_by_names(\@) {
    my $self = shift;
    my @names = @_;
    my @list;
    foreach my $name (@names) {
        my $locus = $self->get_locus_by_id($name);
        push(@list, $locus->get_line());
    }
    return @list;
}

sub write_file_from_loci_ids($\@) {
    my $self = shift;
    my $file_name = shift;
    my @ids = @_;
    open(LOCUS_FILE, ">$file_name")
        or die("Can't open $file_name for writing.\n");
    print LOCUS_FILE $self->header_line();
    my $previous;
    # Get a list of loci instances
    my @loci = map {$self->get_locus_by_id($_) } @ids;
    # Sort, so the loci in the output file are in the same order as
    # in the source file, i.e. that the loci are sorted by position.
    @loci = sort { $a->position() <=> $b->position() } @loci;
    foreach my $locus (@loci) {
        print LOCUS_FILE $locus->get_line($previous) . "\n";
        $previous = $locus;
    }
    close(LOCUS_FILE);
}

sub get_locus_by_id ($) {
    my $self = shift;
    my $id = shift;
    if (not exists $self->{LOCI}{$id}) {
        die("Locus $id not found.\n");
    }
    return $self->{LOCI}{$id};
}

sub get_locus_line_by_id($$) {
    my $self = shift;
    my $id = shift;
    my $previous = shift;
    my $locus = $self->get_locus_by_id($id);
    return $locus->get_line($previous);
}

sub header_line {
    my $self = shift;
    # Do not allow to set the header line as it should only be set in
    # the constructor.
    # if (@_) { $self->{HEADER_LINE} = shift; }
    return $self->{HEADER_LINE};
}

sub get_all_numbers {
    my $self = shift;
    # <=>, the startship operator, ensures numerical sorting.
    return(sort {$a <=> $b} keys %{$self->{NUMBERS}});
}

sub locus_cmp_by_keys {
    my $self = shift;
    my $left = $self->get_locus_by_id($a)->position();
    my $right = $self->get_locus_by_id($b)->position();
    if ($left < $right) {
        return -1;
    } elsif ($left == $right) {
        return 0;
    } else {
        return 1;
    }
}

sub get_all_ids {
    my $self = shift;
    # Ensure they are sorted as in the original file, i.e. according to
    # their numbers.
    # print "Getting all the loci as a sorted list.\n";
    # return sort {$self->locus_cmp_by_keys} keys %{$self->{LOCI}};
    my @numbers = $self->get_all_numbers();
    return map { $self->get_locus_by_number($_)->snp_id() } @numbers;
}

sub get_locus_by_number($) {
    my $self = shift;
    my $number = shift;
    if (not exists $self->{NUMBERS}{$number}) {
        die("Locus number $number not found.\n");
    }
    return $self->{NUMBERS}{$number};
}

sub get_locus_by_position($) {
    my $self = shift;
    my $number = shift;
    if (not exists $self->{POSITIONS}{$number}) {
        die("Locus position '$number' not found.\n");
    }
    return $self->{POSITIONS}{$number};
}

# Get loci names.
# NB: they are returned in the same order they were in the original
# file.
sub get_loci_names {
    my $self = shift;
    my @numbers = $self->get_all_numbers();
    my @names = map { $self->get_locus_by_number($_)->snp_id() } @numbers;
    return @names;
}

# Takes a list of loci IDs as an argument and returns another list of
# IDs, including offset loci.
sub offset($\@) {
    # print "Hapmix::Loci offset()\n";
    my $self = shift;
    my $offset = shift;
    # Initial don't have to be in any order.
    my @initial = @_;
    # Potentials need to be sorted according to their positions,
    # because shortcuts will rely on the ordering.
    my @id_list = $self->get_all_ids();
    my %potential;
    # print "Building potentials()\n";
    foreach my $id (@id_list) {
        $potential{$id} = 1;
    }
    # Target list is initally empty. It is a hash map so loci don't
    # repeat.
    my %offset_ids;
    # For each initial, find its neighbours.
    my @hashes;
    # This loop takes a relatively long time. Perhaps it should display
    # its progress.
    for my $idx (@initial) {
        $self->find_neighbours($offset, $idx, \%potential, \%offset_ids);
    }
    return keys %offset_ids;
}

# Finds neighbours from the potential offset loci.
sub find_neighbours {
    my $self = shift;
    my $offset = shift;
    my $idx = shift;
    my $p_ref = shift;
    my $o_ref = shift;
    # print "find_neighbours($offset, $idx)\n";
    $self->find_neighbours_directioned($offset, $idx, $p_ref, $o_ref, "forward");
    $self->find_neighbours_directioned($offset, $idx, $p_ref, $o_ref, "backward");
}

# Find neighbours by looks at only one direction. It takes advantage
# from the loci positions and numbers to be sorted accordingly.
sub find_neighbours_directioned {
    my $self = shift;
    my $offset = shift;
    my $idx = shift;
    my $p_ref = shift;
    my $o_ref = shift;
    my $direction = shift;
    my $step = 0;

    if ($direction eq "forward") {
        $step = 1;
    } elsif ($direction eq "backward") {
        $step = -1;
    } else {
        die("Direction must be either backward or forward.\n");
    }

    # print "find_neighbours_directioned($offset, $idx, \"$direction\")\n";

    my %potential = %{$p_ref}; # gotcha: this _copies_ the data!
                               # Thou shalt not write to this hash.
                               # It won't be visible outside the
                               # function.
    my $occupant = $self->get_locus_by_id($idx);
    # Get occupant's number.
    my $occupant_no = $occupant->number();
    # Go from this number back, to the lower numbers. Check them all for
    # being in the range. Once a out-of-range locus is spotted, stop.
    my $spotted_outsider = 0;
    my $investigated_locus_no = $occupant_no + $step;
    my $highest_number = $self->{HIGHEST_NUMBER};
    while (
            $investigated_locus_no >= 0
            and $investigated_locus_no <= $highest_number
            and not $spotted_outsider) {
        # Get the locus for investigation.
        my $inv_locus = $self->get_locus_by_number($investigated_locus_no);
        # Check if this locus is in the potential loci set.
        if (exists $potential{$inv_locus->snp_id()}) {
            if ($occupant->within_range($offset, $inv_locus)) {
                # Add to offset loci, delete from potentials, go on.
                $o_ref->{$inv_locus->snp_id()} = 1;
                delete $p_ref->{$inv_locus->snp_id()};
            } else {
                # Outside the range, stop processing.
                $spotted_outsider = 1;
            }
        }
        # If not a potential locus, just go on.
        $investigated_locus_no += $step;
    }
    # print "find_neighbours_directioned(): Finished with ";
    # print $occupant_no - $investigated_locus_no;
    # print " bp offset.\n";
}

sub get_min_max_loci(@) {
    my $self = shift;
    my @ids = @_;
    my @positions = map { $self->get_locus_by_id($_)->position() } @ids;
    my $pos_min = min(@positions);
    my $pos_max = max(@positions);
    my $locus_min = $self->get_locus_by_position($pos_min);
    my $locus_max = $self->get_locus_by_position($pos_max);
    return ($locus_min, $locus_max);
}

sub get_range {
    my $self = shift;
    my $locus_min = shift;
    my $locus_max = shift;
    # Find all the loci that are within the range
    my $number_min = $locus_min->number();
    my $number_max = $locus_max->number();
    # print "no min: '$number_min', no max: '$number_max'\n";
    my @numbers = ($number_min .. $number_max);
    # print "number range: @numbers\n";
    # return IDs by numbers.
    return map { $self->get_locus_by_number($_)->snp_id() } @numbers;
}

sub get_no_loci {
    my $self = shift;
    return $self->{HIGHEST_NUMBER} + 1;
}

# The following line is necessary.
1;
