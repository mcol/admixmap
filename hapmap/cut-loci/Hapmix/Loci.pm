package Hapmix::Loci;
use strict;
sub new {
    my $class = shift;
    my $self = {};
    # Read the data.
    my $locus_file = shift;

    open(LOCI, "<$locus_file")
        or die("Can't open the locus file: $locus_file.\n");
    my @lines = <LOCI>;
    close(LOCI);

    $self->{HEADER_LINE} = shift(@lines);
    # my @headers = split(/\s+/, $self->{HEADER_LINE});

    $self->{LOCI_ARRAY} = \@lines;
    $self->{LOCI} = {};
    foreach my $line (@lines) {
        my @cols = split(" ", $line);
        $self->{LOCI}{$cols[0]} = \@cols;
    }

    bless($self, $class);
    # print "Loci constructed from $locus_file.\n";
    return $self;
}

sub get_lines_by_names(\@) {
    my $self = shift;
    my @names = @_;
    my @list;
    foreach my $name (@names) {
        # print "Locus name: $name\n";
        push(@list, $self->get_locus_as_string($name));
    }
    return @list;
}

sub get_locus_as_list($) {
    my $self = shift;
    my $name = shift;
    return @{$self->{LOCI}{$name}};
}

sub get_locus_as_string($) {
    my $self = shift;
    my $name = shift;
    return join(" ", $self->get_locus_as_list($name));
}

# The following line is necessary.
1;
