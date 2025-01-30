package hit;
require Exporter;

use lib ("$ENV{GEMC}/io");
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(init_hit print_hit);


# Initialize hash maps
sub init_hit
{
	my %hit = ();
	
	# The default value for identifier is "id"
	$hit{"identifiers"} = "id";
	return %hit;
}


# Print hit to TEXT file or upload it onto the DB
sub print_hit
{
	my %configuration = %{+shift};
	my %hit           = %{+shift};
	
	my $varia = $configuration{"variation"};
	my $runno = $configuration{"run_number"};

	# converting the hash maps in local variables
	# (this is necessary to parse the MYSQL command)
	
	my $lname             = trim($hit{"name"});
	my $ldescription      = trim($hit{"description"});
	my $lidentifiers      = trim($hit{"identifiers"});
	my $lSignalThreshold  = trim($hit{"signalThreshold"});
	my $lTimeWindow       = trim($hit{"timeWindow"});
	my $lProdThreshold    = trim($hit{"prodThreshold"});
	my $lMaxStep          = trim($hit{"maxStep"});
	my $lriseTime         = trim($hit{"riseTime"});
	my $lfallTime         = trim($hit{"fallTime"});
	my $lmvToMeV          = trim($hit{"mvToMeV"});
	my $lpedestal         = trim($hit{"pedestal"});
	my $ldelay            = trim($hit{"delay"});

    # after perl 5.10 once can use "state" to use a static variable`
	state $counter_text = 0;
	state $counter_mysql = 0;
	state $counter_sqlite = 0;

	# TEXT Factory
	if($configuration{"factory"} eq "TEXT") {
		my $file = $configuration{"detector_name"}."__hit_".$varia.".txt";
		if($counter_text == 0 || $this_variation ne  $varia)
		{
			`rm -f $file`;
			print "Overwriting if existing: ",  $file, "\n";
			$counter_text = 1;
			$this_variation = $varia;
		}
		
		open(my $info, ">>$file");
		printf $info ("%20s  |", $lname);
		printf $info ("%30s  |", $ldescription);
		printf $info ("%40s  |", $lidentifiers);
		printf $info ("%8s   |", $lSignalThreshold);
		printf $info ("%8s   |", $lTimeWindow);
		printf $info ("%8s   |", $lProdThreshold);
		printf $info ("%8s   |", $lMaxStep);
		printf $info ("%8s   |", $lriseTime);
		printf $info ("%8s   |", $lfallTime);
		printf $info ("%8s   |", $lmvToMeV);
		printf $info ("%8s   |", $lpedestal);
		printf $info ("%8s  \n", $ldelay);
		close($info);
	}
	
	# MYSQL Factory
	if($configuration{"factory"} eq "MYSQL") {

	}

	# SQLITE Factory
	if($configuration{"factory"} eq "SQLITE") {
        my $dbh = open_db(%configuration);
        my $system = $configuration{"detector_name"};

        # first time this module is run, delete everything in geometry table for this variation, system and run number
        if ($counter_sqlite == 0 || $this_variation ne $varia) {
            my $sql = "DELETE FROM hits WHERE system = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all hits for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $mnames_string = "system, variation, run, name, description, identifiers, signalThreshold, timeWindow, prodThreshold, maxStep, riseTime, fallTime, mvToMeV, pedestal, delay";
        # for each name in $mnames_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $mnames_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

                my $sql = "INSERT INTO hits ($mnames_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $lname, $ldescription, $lidentifiers, $lSignalThreshold, $lTimeWindow, $lProdThreshold, $lMaxStep, $lriseTime, $lfallTime, $lmvToMeV, $lpedestal, $ldelay)
        	or die "SQL Error: $DBI::errstr\n";
	}
	
	if($configuration{"verbosity"} > 0)	{
		print "  + Hit $lname uploaded successfully for system: system: $system, variation: \"$varia\", run: $runno \n";
	}
}


1;





