package GXML;

use utils;
use geometry;

sub new
{
	my $class = shift;
	my $self = {
		_dirName => shift,
		_volumes => [],
	};

	bless $self, $class;
	return $self;
}


sub add{
	my $self = shift;
	my %det = %{+shift};

	push(@{$self->{_volumes}}, \%det);

	return %det;
}


sub print {
	my $self = shift;

	open(my $info, ">$self->{_dirName}/cad.gxml");
	printf $info ("<gxml>\n");
	foreach my $det (@{$self->{_volumes}}){
		printf $info ("\t<volume name=\"%s\"", $det->{"name"});
		printf $info (" color=\"%s\"", $det->{"color"});
		printf $info (" material=\"%s\"", $det->{"material"});
		printf $info (" position=\"%s\"", $det->{"pos"});
		printf $info (" rotation=\"%s\"", $det->{"rotation"});
		if($det->{"mother"} ne ""){
			printf $info (" mother=\"%s\"", $det->{"mother"});
		}
		if($det->{"sensitivity"} ne "no"){
			printf $info (" sensitivity=\"%s\"", $det->{"sensitivity"});
			if($det->{"hitType"} ne "no"){
				printf $info (" hitType=\"%s\"", $det->{"hitType"});
			}
			printf $info (" identifiers=\"%s\"", $det->{"identifiers"});
		}
		printf $info (" />\n");
	}
	printf $info ("</gxml>\n");
	close($info);
}


1;
