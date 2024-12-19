package materials;
require Exporter;

use lib ("$ENV{GEMC}/io");
use warnings;
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(init_mat print_mat);


# Initialize hash maps
sub init_mat {
    my %mat = ();

    # The default value for identifier is "id"
    $mat{"description"} = "id";

    # The optical properties are defaulted to none
    # User can define a optical property with arrays of:
    #
    # - photon wavelength (mandatory)
    # - At least one of the following quantities arrays
    $mat{"photonEnergy"} = "none";
    $mat{"indexOfRefraction"} = "none";
    $mat{"absorptionLength"} = "none";
    $mat{"reflectivity"} = "none";
    $mat{"efficiency"} = "none";

    # scintillation specific
    $mat{"fastcomponent"} = "none";
    $mat{"slowcomponent"} = "none";
    $mat{"scintillationyield"} = "-1";
    $mat{"resolutionscale"} = "-1";
    $mat{"fasttimeconstant"} = "-1";
    $mat{"slowtimeconstant"} = "-1";
    $mat{"yieldratio"} = "-1";
    $mat{"rayleigh"} = "none";
    $mat{"birkConstant"} = "-1";

    # mie scattering
    $mat{"mie"} = "none";
    $mat{"mieforward"} = "-1";
    $mat{"miebackward"} = "-1";
    $mat{"mieratio"} = "-1";

    return %mat;
}


# Print material to TEXT file or upload it onto the DB
sub print_mat {
    my %configuration = %{+shift};
    my %mats = %{+shift};

    my $system = $configuration{"detector_name"};
    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    # converting the hash maps in local variables
    # (this is necessary to parse the MYSQL command)
    my $lname = trim($mats{"name"});
    my $ldesc = trim($mats{"description"});
    my $ldensity = trim($mats{"density"});
    my $lncomponents = trim($mats{"ncomponents"});
    my $lcomponents = trim($mats{"components"});

    # optical properties
    my $lphotonEnergy = trim($mats{"photonEnergy"});
    my $lindexOfRefraction = trim($mats{"indexOfRefraction"});
    my $labsorptionLength = trim($mats{"absorptionLength"});
    my $lreflectivity = trim($mats{"reflectivity"});
    my $lefficiency = trim($mats{"efficiency"});

    # scintillation specific
    my $lfastcomponent = trim($mats{"fastcomponent"});
    my $lslowcomponent = trim($mats{"slowcomponent"});
    my $lscintillationyield = trim($mats{"scintillationyield"});
    my $lresolutionscale = trim($mats{"resolutionscale"});
    my $lfasttimeconstant = trim($mats{"fasttimeconstant"});
    my $lslowtimeconstant = trim($mats{"slowtimeconstant"});
    my $lyieldratio = trim($mats{"yieldratio"});
    my $lrayleigh = trim($mats{"rayleigh"});
    my $lbirkConstant = trim($mats{"birkConstant"});

    # mie scattering
    my $lmie = trim($mats{"mie"});
    my $lmieforward = trim($mats{"mieforward"});
    my $lmiebackward = trim($mats{"miebackward"});
    my $lmieratio = trim($mats{"mieratio"});

    # after perl 5.10 once can use "state" to use a static variable`
    state $counter_text = 0;
    state $counter_mysql = 0;
    state $counter_sqlite = 0;
    state $this_variation = "";

    # TEXT Factory
    if ($configuration{"factory"} eq "TEXT") {
        my $file = $configuration{"detector_name"} . "__materials_" . $varia . ".txt";

        if ($counter_text == 0 || $this_variation ne $varia) {
            `rm -f $file`;
            print "Overwriting if existing: ", $file, "\n";
            $counter_text = 1;
            $this_variation = $varia;
        }

        open(my $info, ">>$file");
        printf $info("%20s  |", $lname);
        printf $info("%30s  |", $ldesc);
        printf $info("%10s  |", $ldensity);
        printf $info("%10s  |", $lncomponents);
        printf $info("%50s  |", $lcomponents);

        if ($lphotonEnergy eq "none") {
            printf $info("%5s  |", $lphotonEnergy);
            printf $info("%5s  |", $lindexOfRefraction);
            printf $info("%5s  |", $labsorptionLength);
            printf $info("%5s  |", $lreflectivity);
            printf $info("%5s  |", $lefficiency);

            # scintillation
            printf $info("%5s  |", $lfastcomponent);
            printf $info("%5s  |", $lslowcomponent);
            printf $info("%5s  |", $lscintillationyield);
            printf $info("%5s  |", $lresolutionscale);
            printf $info("%5s  |", $lfasttimeconstant);
            printf $info("%5s  |", $lslowtimeconstant);
            printf $info("%5s  |", $lyieldratio);
            printf $info("%5s  |", $lrayleigh);
            printf $info("%5s  |", $lbirkConstant);
            printf $info("%5s  |", $lmie);
            printf $info("%5s  |", $lmieforward);
            printf $info("%5s  |", $lmiebackward);
            printf $info("%5s  \n", $lmieratio);

        }
        else {
            printf $info("%s  |", $lphotonEnergy);

            # index of refraction
            if ($lindexOfRefraction eq "none") {printf $info("%5s |", $lindexOfRefraction);}
            else {printf $info("%s  |", $lindexOfRefraction);}
            # absorption length
            if ($labsorptionLength eq "none") {printf $info("%5s |", $labsorptionLength);}
            else {printf $info("%s  |", $labsorptionLength);}
            # reflectivity
            if ($lreflectivity eq "none") {printf $info("%5s |", $lreflectivity);}
            else {printf $info("%s  |", $lreflectivity);}
            # efficiency
            if ($lefficiency eq "none") {printf $info("%5s |", $lefficiency);}
            else {printf $info("%s  |", $lefficiency);}

            # scintillation

            # fast component (as function of wavelength)
            if ($lfastcomponent eq "none") {printf $info("%5s |", $lfastcomponent);}
            else {printf $info("%s  |", $lfastcomponent);}
            # slow component (as function of wavelength)
            if ($lslowcomponent eq "none") {printf $info("%5s |", $lslowcomponent);}
            else {printf $info("%s  |", $lslowcomponent);}
            # scintillation yield (constant)
            if ($lscintillationyield eq "-1") {printf $info("%5s |", $lscintillationyield);}
            else {printf $info("%s  |", $lscintillationyield);}
            # resolution scale (constant)
            if ($lresolutionscale eq "-1") {printf $info("%5s |", $lresolutionscale);}
            else {printf $info("%s  |", $lresolutionscale);}
            # fast time (constant)
            if ($lfasttimeconstant eq "-1") {printf $info("%5s |", $lfasttimeconstant);}
            else {printf $info("%s  |", $lfasttimeconstant);}
            # slow time (constant)
            if ($lslowtimeconstant eq "-1") {printf $info("%5s |", $lslowtimeconstant);}
            else {printf $info("%s  |", $lslowtimeconstant);}
            # ratio of yield to total yield for slow component (constant)
            if ($lyieldratio eq "-1") {printf $info("%5s |", $lyieldratio);}
            else {printf $info("%s  |", $lyieldratio);}
            # rayleigh scattering
            if ($lrayleigh eq "none") {printf $info("%5s |", $lrayleigh);}
            else {printf $info("%s  |", $lrayleigh);}
            # Birk constant
            if ($lbirkConstant eq "-1") {printf $info("%5s |", $lbirkConstant);}
            else {printf $info("%s  |", $lbirkConstant);}
            # Mie scattering
            if ($lmie eq "none") {printf $info("%5s |", $lmie);}
            else {printf $info("%s  |", $lmie);}
            if ($lmieforward eq "-1") {printf $info("%5s |", $lmieforward);}
            else {printf $info("%s  |", $lmieforward);}
            if ($lmiebackward eq "-1") {printf $info("%5s |", $lmiebackward);}
            else {printf $info("%s  |", $lmiebackward);}
            if ($lmieratio eq "-1") {printf $info("%5s\n", $lmieratio);}
            else {printf $info("%s \n", $lmieratio);}

        }

        close($info);
    }

    # MYSQL Factory
    if ($configuration{"factory"} eq "MYSQL") {

    }

    # SQLITE Factory
    if ($configuration{"factory"} eq "SQLITE") {
        my $dbh = open_db(%configuration);

        # first time this module is run, delete everything in geometry table for this variation, system and run number
        if ($counter_sqlite == 0 || $this_variation ne $varia) {
            my $sql = "DELETE FROM materials WHERE system = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all materials for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $mnames_string = "system, variation, run, name, description, density, ncomponents, components, photonEnergy, indexOfRefraction, absorptionLength, reflectivity, efficiency, fastcomponent, slowcomponent, scintillationyield, resolutionscale, fasttimeconstant, slowtimeconstant, yieldratio, rayleigh, birkConstant, mie, mieforward, miebackward, mieratio ";

        # for each name in $mnames_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $mnames_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

        my $sql = "INSERT INTO materials ($mnames_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $lname, $ldesc, $ldensity, $lncomponents, $lcomponents, $lphotonEnergy, $lindexOfRefraction, $labsorptionLength, $lreflectivity, $lefficiency, $lfastcomponent, $lslowcomponent, $lscintillationyield, $lresolutionscale, $lfasttimeconstant, $lslowtimeconstant, $lyieldratio, $lrayleigh, $lbirkConstant, $lmie, $lmieforward, $lmiebackward, $lmieratio)
            or die "SQL Error: $DBI::errstr\n";

    }

    if ($configuration{"verbosity"} > 0) {
        print "  + Material $lname uploaded successfully for system: $system, variation: \"$varia\", run: $runno \n";
    }
}

1;
