#!/usr/bin/perl -w

use strict;
use lib ("$ENV{GEMC}/io");
use utils;
use parameters;

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   upload_parameters <config filename> \n";
 	print "   A file named <paramenters.txt> must be present. \n\n";
}

# Make sure the argument list is correct
if( scalar @ARGV != 1) 
{
	help();
	exit;
}

my $config_file = $ARGV[0];
my %configuration = load_configuration($config_file); 

upload_parameters(%configuration);



