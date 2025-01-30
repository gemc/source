package geometry;
require Exporter;

use lib ("$ENV{GEMC}/io");
use warnings;
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(print_det init_det);


# Initialize hash maps
sub init_det {
    my %detector = ();

    # These default value can be left off on the API
    $detector{"description"} = "no description";
    $detector{"pos"} = "0 0 0";
    $detector{"rotation"} = "0 0 0";
    $detector{"dimensions"} = "0";
    $detector{"color"} = "999999";
    $detector{"mfield"} = "no";
    $detector{"ncopy"} = "1";
    $detector{"pMany"} = 1;
    $detector{"exist"} = 1;
    $detector{"visible"} = 1;
    $detector{"style"} = 0;
    $detector{"sensitivity"} = "no";
    $detector{"hit_type"} = "no";
    $detector{"identifiers"} = "no";

    return %detector;
}


# Print detector to TEXT file or upload it onto the DB
sub print_det {
    my %configuration = %{+shift};
    my %det = %{+shift};

    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    # converting the hash maps in local variables
    # (this is necessary to parse the MYSQL command)

    my $lname = trim($det{"name"});
    my $lmother = trim($det{"mother"});
    my $ldescription = trim($det{"description"});
    my $lpos = trim($det{"pos"});
    my $lrotation = trim($det{"rotation"});
    my $lcolor = trim($det{"color"});
    my $ltype = trim($det{"type"});
    my $ldimensions = trim($det{"dimensions"});
    my $lmaterial = trim($det{"material"});
    my $lmfield = trim($det{"mfield"});
    my $lncopy = trim($det{"ncopy"});
    my $lpMany = trim($det{"pMany"});
    my $lexist = trim($det{"exist"});
    my $lvisible = trim($det{"visible"});
    my $lstyle = trim($det{"style"});
    my $lsensitivity = trim($det{"sensitivity"});
    my $lhit_type = trim($det{"hit_type"});
    my $lidentifiers = trim($det{"identifiers"});

    # after 5.10 once can use "state" to use a static variable`
    state $counter_text = 0;
    state $counter_mysql = 0;
    state $counter_sqlite = 0;
    state $this_variation = "";

    # TEXT Factory
    if ($configuration{"factory"} eq "TEXT") {

        my $file = $configuration{"detector_name"} . "__geometry_" . $varia . ".txt";
        if ($counter_text == 0 || $this_variation ne $varia) {
            `rm -f $file`;
            print "Overwriting if existing: ", $file, "\n";
            $counter_text = 1;
            $this_variation = $varia;
        }

        open(my $info, ">>", $file) or die "Could not open file '$file': $!";
        # notice $info( will not work, there has to be a space after $info
        printf $info("%20s  |", $lname);
        printf $info("%20s  |", $lmother);
        printf $info("%30s  |", $ldescription);
        printf $info("%50s  |", $lpos);
        printf $info("%40s  |", $lrotation);
        printf $info("%7s   |", $lcolor);
        printf $info("%20s  |", $ltype);
        printf $info("%60s  |", $ldimensions);
        printf $info("%20s  |", $lmaterial);
        printf $info("%20s  |", $lmfield);
        printf $info("%6s   |", $lncopy);
        printf $info("%6s   |", $lpMany);
        printf $info("%6s   |", $lexist);
        printf $info("%4s   |", $lvisible);
        printf $info("%4s   |", $lstyle);
        printf $info("%20s  |", $lsensitivity);
        printf $info("%20s  |", $lhit_type);
        printf $info("%40s \n", $lidentifiers);
        close($info);
    }

    # MYSQL Factory
    if ($configuration{"factory"} eq "MYSQL") {
    }

    # SQLITEs Factory
    if ($configuration{"factory"} eq "SQLITE") {

        my $dbh = open_db(%configuration);
        my $system = $configuration{"detector_name"};

        # first time this module is run, delete everything in geometry table for this variation, system and run number
        if ($counter_sqlite == 0 || $this_variation ne $varia) {
            my $sql = "DELETE FROM geometry WHERE system = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all volumes for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $names_string = "system, variation, run, name, mother, description, pos, rot, col, type, dimensions, material, magfield, ncopy, pMany, exist, visible, style, sensitivity, hitType, identity";

        # for each name in $names_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $names_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

        my $sql = "INSERT INTO geometry ($names_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $lname, $lmother, $ldescription, $lpos, $lrotation, $lcolor, $ltype, $ldimensions, $lmaterial, $lmfield, $lncopy, $lpMany, $lexist, $lvisible, $lstyle, $lsensitivity, $lhit_type, $lidentifiers)
            or die "Can't execute insert statement: $DBI::errstr";
    }

    if ($configuration{"verbosity"} > 0) {
        print "  + CAD $lname uploaded successfully for system: $system, variation: \"$varia\", run: $runno \n";
    }
}

1;
