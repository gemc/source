package parameters;
require Exporter;

use lib ("$ENV{GEMC}/io");
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(upload_parameter upload_parameters get_parameters get_volumes);


# Utility to upload the parameters from a file
sub upload_parameters {
    my %configuration = %{+shift};

    my $filename = shift;
    my $system = shift;
    my $variation = shift;
    my $runno = shift;
    my $dbhost = $configuration{"dbhost"};
    my $factory = $configuration{"factory"};

    # print configuration
    print "  + Uploading parameters from file: $filename \n";
    print "  + System: $system \n";
    print "  + Variation: $variation \n";
    print "  + Run number: $runno \n";
    print "  + DBhost: $dbhost \n";
    print "  + Factory: $factory \n";


    # Insert parameters into table for the variation, run min, run max. Increment id for this set.
    my $dbh = open_db(%configuration);

    if ($factory eq "MYSQL") {
        my $sql = "LOAD DATA LOCAL INFILE '$filename' \
        INTO TABLE parameters                               \
        FIELDS TERMINATED BY '|'                        \
        LINES TERMINATED BY '\n'
        (parameter_name, value, unit, description, authors, emails, document, variable_name, drawing_author, document_date)
        set variation = '$variation', run = $runno,system = '$system'";
        $dbh->do($sql) or die "Could not execute LOAD DATA INFILE: $DBI::errstr";
    }
    elsif ($factory eq "SQLITE") {

        # first delete all parameters for this system, variation and run number
        my $sql = "DELETE FROM parameters WHERE system = ? and variation = ? and run = ?";
        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $variation, $runno);

        open my $fh, '<', $filename or die "Could not open file '$filename' $!";
        while (my $line = <$fh>) {
            chomp $line;
            my @values = split(/\|/, $line);
            my $pname = trim($values[0]);
            my $pval = $values[1];
            my $puni = $values[2];
            my $pdes = $values[3];
            my $pauts = $values[4];
            my $pema = $values[5];
            my $pdocs = $values[6];
            my $vname = $values[7];
            my $dauth = $values[8];
            my $ddate = $values[9];

            my $mnames_string = "system, variation, run, name, value, vunit, description, authors, emails, document, var_name, doc_author, doc_date";

            # for each name in $mnames_string, we need to add a ? to the values string
            my $qvalues_string = "";
            my @names = split(/\s+/, $mnames_string);
            foreach my $name (@names) {
                $qvalues_string = $qvalues_string . "?, ";
            }
            # remove last comma from $qvalues_string
            $qvalues_string = substr($qvalues_string, 0, -2);

            my $sql = "INSERT INTO parameters ($mnames_string) VALUES ($qvalues_string)";

            my $sth = $dbh->prepare($sql);
            $sth->execute($system, $variation, $runno, $pname, $pval, $puni, $pdes, $pauts, $pema, $pdocs, $vname, $dauth, $ddate)
                or die "Can't execute insert statement: $DBI::errstr";
        }
    }

    $dbh->disconnect();

}

sub upload_parameter {

    if (@_ != 11) {
        print " ERROR: To define a parameter 10 arguments should be passed to <upload_parameter> \n";
    }

    my %configuration = %{+shift};

    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    my $pname = shift; # parameter name
    my $pval = shift;  # parameter value
    my $puni = shift;  # parameter unit
    my $pdes = shift;  # description
    my $pauts = shift; # authors
    my $pema = shift;  # emails
    my $pdocs = shift; # document
    my $vname = shift; # variable name
    my $dauth = shift; # drawing author
    my $ddate = shift; # document date

    # after perl 5.10 once can use "state" to use a static variable`
    state $counter_text = 0;
    state $counter_mysql = 0;
    state $counter_sqlite = 0;
    state $this_variation = "";

    # TEXT Factory
    if ($configuration{"factory"} eq "TEXT" || $this_variation ne $varia) {
        my $file = $configuration{"detector_name"} . "__parameters_" . $varia . ".txt";
        if ($counter_text == 0) {
            `rm -f $file`;
            print "Overwriting if existing: ", $file, "\n";
            $counter_text = 1;
            $this_variation = $varia;
        }

        open(my $info, ">>$file");
        printf $info("%20s  |", $pname);
        printf $info("%20s  |", $pval);
        printf $info("%20s  |", $puni);
        printf $info("%50s  |", $pdes);
        printf $info("%20s  |", $pauts);
        printf $info("%20s  |", $pema);
        printf $info("%20s  |", $pdocs);
        printf $info("%20s  |", $vname);
        printf $info("%20s  |", $dauth);
        printf $info("%20s  \n", $ddate);
        close($info);
    }

    # MYSQL Factory
    if ($configuration{"factory"} eq "MYSQL") {

    }

    # SQLITE Factory
    if ($configuration{"factory"} eq "SQLITE") {
        my $dbh = open_db(%configuration);
        my $system = $configuration{"detector_name"};

        # first time this module is run, delete everything in geometry table for this variation, system and run number
        if ($counter_sqlite == 0 || $this_variation ne $varia) {
            my $sql = "DELETE FROM parameters WHERE system = ? and variation = ? and run = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all parameters for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $mnames_string = "system, variation, run, parameter_name, value, unit, description, authors, emails, document, variable_name, drawing_author, document_date";

        # for each name in $mnames_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $mnames_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

        my $sql = "INSERT INTO parameters ($mnames_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $pname, $pval, $puni, $pdes, $pauts, $pema, $pdocs, $vname, $dauth, $ddate)
            or die "Can't execute insert statement: $DBI::errstr";
    }

    if ($configuration{"verbosity"} > 0) {
        print "  + variable $lname uploaded successfully for variation \"$varia\" \n";
    }

}


# Utility to get a hash map with the parameters
sub get_parameters {
    my (%configuration) = @_;
    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    my %parameters = ();

    # Text Factory. The parameter file is assumed to be present
    # and named "parameters.txt"
    if ($configuration{"factory"} eq "TEXT") {
        my $file = $configuration{"detector_name"} . "__parameters_" . $varia . ".txt";
        open(FILE, $file) or die("Open failed on file $file: $!");
        my @lines = <FILE>;
        close(FILE);
        foreach my $line (@lines) {
            my @numbers = split(/[|]+/, $line);
            my ($pnam, $pval) = @numbers;
            $parameters{trim($pnam)} = trim($pval);
        }
    }

    if ($configuration{"factory"} eq "MYSQL") {
    }

    if ($configuration{"factory"} eq "SQLITE") {
        my $dbh = open_db(%configuration);

        # Get the correct run ranges and the latest version
        my $query = $dbh->prepare("SELECT * from parameters where variation = '$varia' and run = '$runno' ");
        $query->execute();

        while (my @data = $query->fetchrow_array()) {
            my $pnam = $data[4];
            my $pval = $data[5];
            $parameters{trim($pnam)} = trim($pval);
        }

        $dbh->disconnect();
    }

    if ($configuration{"verbosity"} > 0) {
        foreach my $key (keys %parameters) {
            print " * Parameter \"$key\" loaded with value: $parameters{$key} \n";
        }
    }

    print "\n";
    return %parameters;
}

# Utility to get a hash maps with the volumes
# Subroutine to read txt file with volumes from COATJAVA FTOF factory
sub get_volumes {
    my (%configuration) = @_;
    my $varia = $configuration{"variation"};

    # Hash maps to populate from volumes files
    my %mothers = ();
    my %positions = ();
    my %rotations = ();
    my %types = ();
    my %dimensions = ();
    my %ids = ();

    # Text Factory. The volumes file is assumed to be present
    # and named "volumes.txt"
    # If it is not present run COATJAVA Detector Factory to create it
    my $file = $configuration{"detector_name"} . "__volumes_" . $varia . ".txt";
    open(FILE, $file) or die("Open failed on file $file: $! (Run factory.groovy to create this file)");
    my @lines = <FILE>;
    close(FILE);
    foreach my $line (@lines) {
        my @vvalues = split('[|]+', $line);

        # Assign fields to corresponding hash maps
        $vnam = trim($vvalues[0]);
        $mothers{$vnam} = trim($vvalues[1]);
        $positions{$vnam} = trim($vvalues[2]);
        $rotations{$vnam} = trim($vvalues[3]);
        $types{$vnam} = trim($vvalues[4]);
        $dimensions{$vnam} = trim($vvalues[5]);
        $ids{$vnam} = trim($vvalues[6]);
    }

    if ($configuration{"verbosity"} > 1) {
        foreach my $key (keys %dimensions) {
            print " * Parameter \"$key\" loaded with mother: $mothers{$key} \n";
            print " * Parameter \"$key\" loaded with position: $positions{$key} \n";
            print " * Parameter \"$key\" loaded with rotation: $rotations{$key} \n";
            print " * Parameter \"$key\" loaded with type: $types{$key} \n";
            print " * Parameter \"$key\" loaded with dimension: $dimensions{$key} \n";
        }
    }

    print "\n";
    return (\%mothers, \%positions, \%rotations, \%types, \%dimensions, \%ids);
}

