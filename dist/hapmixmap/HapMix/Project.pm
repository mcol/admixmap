#!/usr/bin/perl
package HapMix::Project;

use strict;

use Time::Local;
use HapMix::Common;
use Config::IniFiles;
use File::Copy;

sub new {
    my $class = shift;
    my $self = {};
    $self->{NAME} = shift;
    my $config_file = "$self->{NAME}/project.ini";
    if (-r $config_file) {
        $self->{CONFIG} = Config::IniFiles->new(-file => $config_file);
    } else {
        # This project doesn't have a config.
        $self->{CONFIG} = Config::IniFiles->new();
        $self->{CONFIG}->SetFileName($config_file);
    }
    $self->{COMMON_DATA} = HapMix::Common->new();
    bless($self, $class);
    return $self;
}

sub initialize {
    my $self = shift;
    my $population = shift;
    my $project_name = $self->{NAME};
    mkdir($project_name, 0755)
        or die("Can't create directory $project_name");
    my @subdirs = keys %{$self->{COMMON_DATA}->{DIRECTORIES}};
    for my $subdir_key (@subdirs) {
        my $subdir = $self->{COMMON_DATA}->{DIRECTORIES}{$subdir_key};
        mkdir("$project_name/$subdir", 0755)
            or die("Can't create directory $project_name/$subdir.\n");
    }
    $self->{CONFIG}->AddSection("project");
    $self->{CONFIG}->newval("project", "name", $project_name);
    $self->{CONFIG}->newval("project", "population", $population);
    $self->{CONFIG}->newval("project", "created_on", scalar(localtime()));
    $self->{CONFIG}->RewriteConfig();
}

sub get_dir {
    my $self = shift;
    return $self->{CONFIG}->val("project", "name");
}

sub get_source_file_name {
    my $self = shift;
    return $self->get_dir()
        . "/"
        . $self->{COMMON_DATA}->{DIRECTORIES}{source} 
        . "/"
        . $self->{COMMON_DATA}->{SRC_FILE};
}

sub get_source {
    my $self = shift;
    my $source_file = shift;
    my $dst = $self->get_source_file_name;
    copy($source_file, $dst)
        or die("Can't copy $source_file to $dst.\n");
}

1;
