package Hapmix::Masking;
use strict;

# Get two lines of an individual by individual no
sub get_indiv_lines {
	my $indiv_no = shift;
	my $lines = shift;
	my $idx = $indiv_no * 2;
	return ($lines[$idx], $lines[$idx + 1]);
}

# Arrays for individuals
my (@indivs_cc, @indivs_train);
@indivs_cc = @indiv_idx[0 .. ($no_masked - 1)];
@indivs_train = @indiv_idx[$no_masked .. $#indiv_idx];

sub compose_diploid(\@\@$) {
	# my $lines_ref = @_;
	# my @lines = @$lines_ref;
	my ($cols1_ref, $cols2_ref, $limit) = @_;
	my @outlist;
	my @cols1 = @$cols1_ref;
	my @cols2 = @$cols2_ref;
	$outlist[0] = shift(@cols1);
	shift(@cols2);
	foreach my $i (0 .. ($limit - 1)) {
		$outlist[$i + 1] = '"' . $cols1[$i] . "," . $cols2[$i] . '"';
	};
	return join("\t", @outlist) . "\n";
}

sub write_idx_to_file(\$\$$\@\@$) {
	my ($filename_ref, $header_ref, $format, $indices_ref, $wlines_ref, $limit_loci) = @_;
	my $filename = $$filename_ref;
	my $header = $$header_ref;
	my @indices = @$indices_ref;
	my @wlines = @$wlines_ref;
	my $limit = 0;
	# Valid formats
	my %formats = (
		'diploid' => 1,
		'haploid' => 1);
	if (!($formats{$format})) {
		die("Format ($format) should be either diploid or haploid.");
	}
	open(OUT_FILE, ">$filename");
	my @firstline = split(/\s+/, $header);
	if ($format eq "diploid") {
		$firstline[0] = "Individ";
	}
	$limit = $limit_loci ? $limit_loci : ($#firstline - 1);
    # If the limit is less than the number of loci, use the number of
    # loci.
    if ($limit_loci and $limit_loci < $#firstline) {
        $limit =  $limit_loci;
    } else {
        $limit =  ($#firstline - 1);
    }
	print OUT_FILE join("\t", @firstline[0 .. $limit]) . "\n";
	foreach my $idx (@indices) {
		my @ind_lines = get_indiv_lines($idx, @wlines);
		# print substr($ind_lines[0], 0, 72) . "\n";
		# print substr($ind_lines[1], 0, 72) . "\n";
		@cols1 = split(/\s+/, $ind_lines[0]);
		@cols2 = split(/\s+/, $ind_lines[1]);
		# Check for label correctness, every two rows:
		# LABEL
		# LABEL_2
		if (!($cols2[0] eq ($cols1[0] . "_2"))) {
			die("Wrong gamete labels: $cols1[0], $cols2[0].");
		}
		if ($format eq 'haploid') {
			print OUT_FILE join(" ", @cols1[0 .. $limit]) . "\n";
			print OUT_FILE join(" ", @cols2[0 .. $limit]) . "\n";
		} elsif ($format eq 'diploid') {
			print OUT_FILE compose_diploid(@cols1, @cols2, $limit);
			# print compose_diploid(@cols1, @cols2);
		} else {
			die("This shouldn't happen.");
		}
	}
	close(OUT_FILE);
}

sub rewrite_loci($$$) {
	my ($in, $out, $limit) = @_;
	open(IN_LOCUSFILE, "<$in")
		or die ("Can't open the input locus file '$in' for reading.");
	open(OUT_LOCUSFILE, ">$out")
		or die ("Can't open the output locus file '$out' for writing.");
	my $counter = 0;
	while (<IN_LOCUSFILE>) {
		if ($limit && $counter <= $limit) {
			print OUT_LOCUSFILE $_;
		}
		$counter++;
	}
	close(IN_LOCUSFILE);
	close(OUT_LOCUSFILE);
}



1;
