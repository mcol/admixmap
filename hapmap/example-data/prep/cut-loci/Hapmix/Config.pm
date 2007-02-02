package Hapmix::Config;
use strict;

# Configuration for the package. I don't know if it's the best way,
# but...
sub new {
    my $class = shift;
    my $self = {};
    $self->{SEPARATOR} = "\t";
    bless($self, $class);
    return $self;
}

sub separator {
    my $self = shift;
    if (@_) {
        $self->{SEPARATOR} = shift;
    }
    return $self->{SEPARATOR};
}

sub sep {
    my $self = shift;
    return $self->separator(@_);
}


1;
