package math;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw($pi $inches asind phi theta rad deg sq);

our $pi     = 3.141592653589793238;
$inches = 2.54;
use Math::Trig;

sub asind {asin($_[0])*180/$pi}

sub phi
{
	my $x = $_[0];
	my $y = $_[1];
	my $z = $_[2];
	
	my $phi = 0.0;
	if($y > 0.0)
	{
		$phi = acos($x/sqrt($x*$x+$y*$y));
	}
	if($y < 0.0)
	{
    $phi =  2*3.141592654 - acos($x/sqrt($x*$x+$y*$y));
	}
	if($y == 0.0 && $x<0)
	{
    $phi =  2*3.141592654;
	}
	return $phi;
}

sub theta
{
	my $x = $_[0];
	my $y = $_[1];
	my $z = $_[2];
	
	my $theta = 0.0;
	
	if($x*$x+$y*$y+$z*$z)
	{
		$theta =  acos($z/sqrt($x*$x+$y*$y+$z*$z));
	}
	
	return $theta;
}

sub rad
{
	return $_[0]*$pi/180.0;
}

sub deg
{
	return $_[0]*180.0/$pi;
}

sub sq   { $_[0] * $_[0] }


1;





