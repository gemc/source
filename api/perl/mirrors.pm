package mirrors;
require Exporter;

use lib ("$ENV{GEMC}/io");
use utils;
use 5.010;

@ISA = qw(Exporter);
@EXPORT = qw(init_mir print_mir);


# Available finish in materials/include/G4OpticalSurface.hh:
#
# polished,                    // smooth perfectly polished surface
# polishedfrontpainted,        // smooth top-layer (front) paint
# polishedbackpainted,         // same is 'polished' but with a back-paint
#
# ground,                      // rough surface
# groundfrontpainted,          // rough top-layer (front) paint
# groundbackpainted,           // same as 'ground' but with a back-paint
#
# polishedlumirrorair,         // mechanically polished surface, with lumirror
# polishedlumirrorglue,        // mechanically polished surface, with lumirror & meltmount
# polishedair,                 // mechanically polished surface
# polishedteflonair,           // mechanically polished surface, with teflon
# polishedtioair,              // mechanically polished surface, with tio paint
# polishedtyvekair,            // mechanically polished surface, with tyvek
# polishedvm2000air,           // mechanically polished surface, with esr film
# polishedvm2000glue,          // mechanically polished surface, with esr film & meltmount
#
# etchedlumirrorair,           // chemically etched surface, with lumirror
# etchedlumirrorglue,          // chemically etched surface, with lumirror & meltmount
# etchedair,                   // chemically etched surface
# etchedteflonair,             // chemically etched surface, with teflon
# etchedtioair,                // chemically etched surface, with tio paint
# etchedtyvekair,              // chemically etched surface, with tyvek
# etchedvm2000air,             // chemically etched surface, with esr film
# etchedvm2000glue,            // chemically etched surface, with esr film & meltmount
#
# groundlumirrorair,           // rough-cut surface, with lumirror
# groundlumirrorglue,          // rough-cut surface, with lumirror & meltmount
# groundair,                   // rough-cut surface
# groundteflonair,             // rough-cut surface, with teflon
# groundtioair,                // rough-cut surface, with tio paint
# groundtyvekair,              // rough-cut surface, with tyvek
# groundvm2000air,             // rough-cut surface, with esr film
# groundvm2000glue             // rough-cut surface, with esr film & meltmount


# Available models in materials/include/G4OpticalSurface.hh:
#
# glisur,                      // original GEANT3 model
# unified,                     // UNIFIED model
# LUT                          // Look-Up-Table model


# Available surface types in materials/include/G4SurfaceProperty.hh
#
# dielectric_metal,            // dielectric-metal interface
# dielectric_dielectric,       // dielectric-dielectric interface
# dielectric_LUT,              // dielectric-Look-Up-Table interface
# firsov,                      // for Firsov Process
# x_ray                        // for x-ray mirror process


# Border Volume Types:
#
# SkinSurface: surface of a volume
# Border Surface: surface between two volumes (second volume must exist)


# Initialize hash maps
sub init_mir {
    my %mir = ();

    # The default value for identifier is "id"
    $mir{"description"} = "id";

    # The optical properties are defaulted to none
    # User can define a optical property with arrays of:
    #
    # These properties are mandatory
    $mir{"type"} = "mustBeDefined";
    $mir{"finish"} = "mustBeDefined";
    $mir{"model"} = "mustBeDefined";
    $mir{"border"} = "mustBeDefined";

    # If materialOptProperties is defined, use the material
    # optical properties instead of defining new ones
    $mir{"maptOptProps"} = "none";

    # - At least one of the following quantities arrays
    $mir{"photonEnergy"} = "none";
    $mir{"indexOfRefraction"} = "none";
    $mir{"reflectivity"} = "none";
    $mir{"efficiency"} = "none";
    $mir{"specularlobe"} = "none";
    $mir{"specularspike"} = "none";
    $mir{"backscatter"} = "none";
    $mir{"sigmaAlhpa"} = "-1";

    return %mir;
}

# Print mirror to TEXT file or upload it onto the DB
sub print_mir {
    my %configuration = %{+shift};
    my %mirs = %{+shift};

    my $system = $configuration{"detector_name"};
    my $varia = $configuration{"variation"};
    my $runno = $configuration{"run_number"};

    # converting the hash maps in local variables
    # (this is necessary to parse the MYSQL command)
    my $lname = trim($mirs{"name"});
    my $ldesc = trim($mirs{"description"});

    # Mirror mandatory properties
    my $ltype = trim($mirs{"type"});     # see above surface types
    my $lfinish = trim($mirs{"finish"}); # see above finishes
    my $lmodel = trim($mirs{"model"});   # see above model types
    my $lborder = trim($mirs{"border"}); # see above border types

    # If materialOptProperties is defined, use the material
    # optical properties instead of defining new ones
    my $lmat_opt_props = trim($mirs{"maptOptProps"});

    my $lphoton_energy = trim($mirs{"photonEnergy"});
    my $lrefraction_index = trim($mirs{"indexOfRefraction"});
    my $lreflectivity = trim($mirs{"reflectivity"});
    my $lefficiency = trim($mirs{"efficiency"});
    my $specular_lobee = trim($mirs{"specularlobe"});
    my $lspecular_spike = trim($mirs{"specularspike"});
    my $lbackscatter = trim($mirs{"backscatter"});
    my $lsigma_alpha = trim($mirs{"sigmaAlhpa"});

    if ($lmat_opt_props eq "none") {
        if ($lphoton_energy eq "none") {
            print " !! Error: there is no material with optical properties associated with this mirror.\n";
            print " !! Optical properties must be defined.\n";
        }
    }

    # after perl 5.10 once can use "state" to use a static variable`
    state $counter_text = 0;
    state $counter_mysql = 0;
    state $counter_sqlite = 0;
    state $this_variation = "";

    # TEXT Factory
    if ($configuration{"factory"} eq "TEXT") {
        my $file = $configuration{"detector_name"} . "__mirrors_" . $varia . ".txt";

        if ($counter_text == 0 || $this_variation ne $varia) {
            `rm -f $file`;
            print "Overwriting if existing: ", $file, "\n";
            $counter_text = 1;
            $this_variation = $varia;
        }

        open(my $info, ">>$file");
        printf $info ("%20s  |", $lname);
        printf $info ("%30s  |", $ldesc);

        if ($ltype eq "mustBeDefined") {print "Error: type undefined.\n";}
        if ($lfinish eq "mustBeDefined") {print "Error: finish undefined.\n";}
        if ($lmodel eq "mustBeDefined") {print "Error: model undefined.\n";}
        if ($lborder eq "mustBeDefined") {print "Error: border undefined.\n";}

        printf $info ("%24s  |", $ltype);
        printf $info ("%20s  |", $lfinish);
        printf $info ("%10s  |", $lmodel);
        printf $info ("%25s  |", $lborder);
        printf $info ("%25s  |", $lmat_opt_props);


        # photon energy
        if ($lphoton_energy eq "none") {printf $info ("%5s |", $lphoton_energy);}
        else {printf $info ("%s  |", $lphoton_energy);}
        # index of refraction
        if ($lrefraction_index eq "none") {printf $info ("%5s |", $lrefraction_index);}
        else {printf $info ("%s  |", $lrefraction_index);}
        # reflectivity
        if ($lreflectivity eq "none") {printf $info ("%5s |", $lreflectivity);}
        else {printf $info ("%s  |", $lreflectivity);}
        # efficiency
        if ($lefficiency eq "none") {printf $info ("%5s |", $lefficiency);}
        else {printf $info ("%s  |", $lefficiency);}
        # specularlobe
        if ($specular_lobee eq "none") {printf $info ("%5s |", $specular_lobee);}
        else {printf $info ("%s  |", $specular_lobee);}
        # specularspike
        if ($lspecular_spike eq "none") {printf $info ("%5s |", $lspecular_spike);}
        else {printf $info ("%s  |", $lspecular_spike);}
        # backscatter
        if ($lbackscatter eq "none") {printf $info ("%5s |", $lbackscatter);}
        else {printf $info ("%s  |", $lbackscatter);}
        # sigmaAlhpa
        if ($lsigma_alpha eq "-1") {printf $info ("%5s\n", $lsigma_alpha);}
        else {printf $info ("%s \n", $lsigma_alpha);}

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
            my $sql = "DELETE FROM mirrors WHERE system = ?";
            my $sth = $dbh->prepare($sql);
            $sth->execute($system);
            print "   > Deleted all mirrors for system $system \n";
            $counter_sqlite = 1;
            $this_variation = $varia;
        }

        my $mnames_string = "system, variation, run, name, description, type, finish, model, border, mat_opt_props, photon_energy, refraction_index, reflectivity, efficiency, specular_lobe, specular_spike, backscatter, sigma_alpha";

        # for each name in $mnames_string, we need to add a ? to the values string
        my $qvalues_string = "";
        my @names = split(/\s+/, $mnames_string);
        foreach my $name (@names) {
            $qvalues_string = $qvalues_string . "?, ";
        }
        # remove last comma from $qvalues_string
        $qvalues_string = substr($qvalues_string, 0, -2);

        my $sql = "INSERT INTO mirrors ($mnames_string) VALUES ($qvalues_string)";

        my $sth = $dbh->prepare($sql);
        $sth->execute($system, $varia, $runno, $lname, $ldesc, $ltype, $lfinish, $lmodel, $lborder, $lmat_opt_props, $lphoton_energy, $lrefraction_index, $lreflectivity, $lefficiency, $specular_lobee, $lspecular_spike, $lbackscatter, $lsigma_alpha)
            or die "SQL Error: $DBI::errstr\n";
    }

    if ($configuration{"verbosity"} > 0) {
        print "  + Mirror $lname uploaded successfully for system: $system, variation: \"$varia\", run: $runno \n";
    }
}

1;





