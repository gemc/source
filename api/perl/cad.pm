package cad;
require Exporter;

use lib ("$ENV{GEMC}/io");
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(init_cad print_cad);


# Initialize hash maps
sub init_cad {
    my %cad = ();

    # default values
    $cad{"cad_subdir"} = "cad";
    $cad{"sensitivity"} = "no";
    $cad{"hit_type"} = "no";
    $cad{"identifiers"} = "no";
    $cad{"visible"} = 1;
    $cad{"style"} = 1;
    $cad{"position"} = "0*mm, 0*mm, 0*mm";
    $cad{"rotation"} = "0*deg, 0*deg, 0*deg";
    $cad{"mfield"} = "na";
    $cad{"mother"} = "root";
    $cad{"color"} = "778899";
    return %cad;
}


# Print cad to TEXT file or upload it onto the DB
sub print_cad {
    my %configuration = %{+shift};
    my %cad = %{+shift};

    my $system = $configuration{"detector_name"}."_cad";
    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    # converting the hash maps in local variables
    # (this is necessary to parse the MYSQL command)

    my $lname = trim($cad{"name"});
    my $lcad_subdir = trim($cad{"cad_subdir"});
    my $lsensitivity = trim($cad{"sensitivity"});
    my $lhit_type = trim($cad{"hit_type"});
    my $lidentifiers = trim($cad{"identifiers"});
    my $lvisible = trim($cad{"visible"});
    my $lstyle = trim($cad{"style"});
    my $lposition = trim($cad{"position"});
    my $lrotation = trim($cad{"rotation"});
    my $lmfield = trim($cad{"mfield"});
    my $lmother = trim($cad{"mother"});
    my $lmaterial = trim($cad{"material"});
    my $lcolor = trim($cad{"color"});


    # after perl 5.10 once can use "state" to use a static variable`
    state $counter_text = 0;
    state $counter_mysql = 0;
    state $counter_sqlite = 0;
    state $this_variation = "";

    # MYSQL Factory
    if ($configuration{"factory"} eq "MYSQL") {
    }

    # SQLITE Factory
    if ($configuration{"factory"} eq "SQLITE") {

        my $dbh = open_db(%configuration);

        # first time this module is run, delete everything in geometry table for this variation, system and run number
        if ($counter_sqlite == 0 || $this_variation ne $varia) {
            my $sql = "DELETE FROM cad WHERE system = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all cad entries for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $mnames_string = "system, variation, run, name, cad_subdir, sensitivity, hit_type, identifiers, visible, style, position, rotation, mfield, mother, material, color";

        # for each name in $mnames_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $mnames_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

        my $sql = "INSERT INTO cad ($mnames_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $lname, $lcad_subdir, $lsensitivity, $lhit_type, $lidentifiers, $lvisible, $lstyle, $lposition, $lrotation, $lmfield, $lmother, $lmaterial, $lcolor)
            or die "SQL Error: $DBI::errstr\n";
    }

    if ($configuration{"verbosity"} > 0) {
        print "  + CAD $lname uploaded successfully for system: $system, variation: \"$varia\", run: $runno \n";
    }
}

1;





