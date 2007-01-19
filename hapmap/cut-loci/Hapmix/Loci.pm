package Hapmix::Loci;
use strict;
use Hapmix::Locus;

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
    foreach my $idx (@ids) {
        print LOCUS_FILE $self->get_locus_line_by_id($idx, $previous) . "\n";
        $previous = $self->get_locus_by_id($idx);
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
    # <=>, startship operator, ensures numerical sorting.
    return(sort {$a <=> $b} keys %{$self->{NUMBERS}});
}

sub get_locus_by_number($) {
    my $self = shift;
    my $number = shift;
    if (not exists $self->{NUMBERS}{$number}) {
        die("Locus number $number not found.\n");
    }
    return $self->{NUMBERS}{$number};
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

# The following line is necessary.
1;
