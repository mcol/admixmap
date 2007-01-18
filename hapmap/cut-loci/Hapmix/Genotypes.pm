package Hapmix::Genotypes;
use strict;
sub new {
    my $class = shift;
    my $self = {};
    # Read the data.
    my $genotypes_file = shift;

    open(GENOTYPES, "<$genotypes_file")
        or die("Can't open the genotypes file: $genotypes_file.\n");
    # print "Reading the array.\n";
    my @lines = <GENOTYPES>;
    # print "Reading array finished.\n";
    $self->{INDIVS} = {};
    # print "Headers: '$lines[0]'\n";
    my $headers = shift(@lines);
    $self->{HEADER_LINE} = $headers;
    # print "Headers: '$headers'\n";
    my @hdr_arr = split(/\s+/, $headers);
    # first column is "Indiv" anyway
    shift(@hdr_arr);
    $self->{LOCI} = \@hdr_arr;
    # Store indices of loci
    $self->{LOCI_INDICES} = {};
    my $counter = 0;
    foreach my $locus (@hdr_arr) {
        $self->{LOCI_INDICES}{$locus} = $counter++;
    }
    $self->{GENOTYPES_ARRAY} = \@lines;
    # foreach my $line ($self->{GENOTYPES_ARRAY}) {
    $counter = 0;
    foreach my $line (@lines) {
        my @cols = split(/\s+/, $line);
        my $id = shift(@cols);
        undef @cols;
        # $self->{INDIVS}{$id} = \@cols;
        $self->{INDIVS}{$id} = $line;
        $counter++;
        # print "Splitting rows. Count: $counter of " . scalar(@lines) . "\n";
    }
    # print "\n";
    close(GENOTYPES);

    bless($self, $class);
    # print "Genotypes constructed from $genotypes_file.\n";
    return $self;
}

sub get_headers {
    my $self = shift;
    return $self->{HEADER_LINE};
}

sub get_loci {
    my $self = shift;
    return @{$self->{LOCI}};
}

sub get_first_header {
    my $self = shift;
    return @{$self->{HEADER_LINE}}[0];
}

sub get_last_header {
    my $self = shift;
    return @{$self->{HEADER_LINE}}[$#{$self->{HEADER_LINE}}];
}

sub get_column_idx_by_name($) {
    my $self = shift;
    my $name = shift;
    if (not defined $self->{LOCI_INDICES}{$name}) {
        die("Column $name not known.\n");
    }
    return $self->{LOCI_INDICES}{$name};
}

sub get_column_indices_by_names(\@) {
    my $self = shift;
    my @names = @_;
    my @indices;
    # print join(" ", @names) . "\n";
    # die("");
    foreach my $name (@names) {
        my $col_idx = $self->get_column_idx_by_name($name);
        # print "Name: '$name', Col idx: '$col_idx'\n";
        push(@indices, $col_idx);
    }
    # print join(" ", @indices) . "\n";
    # die("");
    return @indices;
}

# Return a list of strings, to be written as a table.
# Needs column indices as parameter.
sub get_columns_by_indices(\@) {
    my $self = shift;
    # my @column_ids = @_;
    # my @indices = $self->get_column_indices_by_names(@column_ids);
    my @indices = @_;
    my @cols;
    # For each individual
    foreach my $key (sort keys(%{$self->{INDIVS}})) {
        # print "Indiv: $key\n";
        # print join(" ", @indices) . "\n";
        if (not exists $self->{INDIVS}{$key}) {
            die("Individual $key not known!\n");
        }
        # print "Indiv: '" . $self->{INDIVS}{$key} . "'";
        my @indiv_line = split(/\s+/, $self->{INDIVS}{$key});
        my $data = join(" ", @indiv_line[@indices]);
        $data = $key . " " . $data;
        push(@cols, $data);
    }
    return @cols;
}

# Writes file with given columns
sub write_cols_to_file($\@) {
    my $self = shift;
    my ($file_name, $indices_ref) = @_;
    # my @cols = @{shift()};
    my @indices = @$indices_ref;
    # print join(" ", @indices) . "\n";
    # die("er.");
    open(OUTFILE, ">$file_name")
        or die("Can't open $file_name for writing.");
    my @head_arr = split(/\s+/, $self->{HEADER_LINE});
    print OUTFILE "Individ " . join(" ", @head_arr[@indices]) . "\n";
    foreach my $line ($self->get_columns_by_indices(@indices)) {
        print OUTFILE $line . "\n";
    }
    close(OUTFILE);
}

sub write_file_from_loci_ids($\@) {
    my $self = shift;
    my $file_name = shift;
    my @ids = @_;
    my @indices = $self->get_column_indices_by_names(@ids);
    # print join(" ", @indices) . "\n";
    # die("er.");
    $self->write_cols_to_file($file_name, \@indices);
}

# The following line is necessary.
1;
