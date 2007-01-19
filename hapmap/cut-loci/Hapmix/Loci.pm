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
    for my $idx (@initial) {
        $self->find_neighbours($offset, $idx, \%potential, \%offset_ids);
    }
    return keys %offset_ids;
}

# Returns a list of references to hashes with updated values.
# FIXME: This function is quite inefficient. It always goes through the
# full list of potential neighbours.
sub find_neighbours {
    my $self = shift;
    my $offset = shift;
    my $idx = shift;
    my $p_ref = shift;
    my $o_ref = shift;
    # print "find_neighbours($offset, $idx)\n";
    my $occupant = $self->get_locus_by_id($idx);
    foreach my $key (keys %{$p_ref}) {
        my $locus = $self->get_locus_by_id($key);
        if ($occupant->within_range($offset, $locus)) {
            $o_ref->{$key} = 1;
            delete $p_ref->{$key};
        }
    }
}

sub get_min_max_loci(@) {
    my $self = shift;
    my @ids = @_;
    # print "ids: @ids\n";
    my @positions = map { $self->get_locus_by_id($_)->position() } @ids;
    # my @positions = keys %{$self->{LOCI}->{POSITIONS}};
    # print "positions: @positions\n";
    my $pos_min = min(@positions);
    my $pos_max = max(@positions);
    # print "min: '$pos_min', max: '$pos_max'\n";
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

# The following line is necessary.
1;
