#! /usr/bin/env perl

use strict;
use warnings;

use Tk;
use Tk::NoteBook;
use Tk::DialogBox;
use Tk::LabEntry;
use XML::Simple;
use Carp;
use Data::Dumper;

use MTToolbox::MTParams 4.1.2;
use MTToolbox::MTEnv 4.1.2;


# Subroutines
sub run;
sub add_sample;
sub remove_sample;
sub add_sample_batch;
sub reset_GUI;
sub guess_file_name;
sub guess_barcode;
sub recalculate_p_nums;
sub _handle_errors;
sub _translate_seq_type;
sub double_click;
sub dirDialog;
sub fileDialog;

sub _output_root_field;
sub _parallelize_field;
sub _fwd_read_file_field;
sub _rev_read_file_field;
sub _fwd_primer_field;
sub _rev_primer_field;
sub _seq_type_field;
sub _keep_tmp_files_field;
sub _fwd_mer_len_field;
sub _rev_mer_len_field;
sub _fwd_max_shifts_field;
sub _rev_max_shifts_field;
sub _fwd_linker_field;
sub _rev_linker_field;
sub _min_con_depth_field;
sub _diginorm_max_field;
sub _con_algo_field;
sub _min_overlap_field;
sub _max_overlap_field;
sub _mismatch_ratio_field;
sub _phred_offset_field;
sub _avg_read_len_field;
sub _avg_frag_len_field;
sub _sd_field;
sub _con_trim_to_base_field;
sub _con_min_len_field;
sub _con_min_avg_qual_field;
sub _con_allow_ambig_field;
sub _con_min_c_score_field;
sub _SRC_trim_to_base_field;
sub _SRC_min_len_field;
sub _SRC_min_avg_qual_field;
sub _SRC_allow_gaps_field;
sub _SRC_allow_ambig_field;




# Some Global Variables
my $list; # a list object that goes on the GUI.
my @sample_arr = ();  # stores all info for each sample in array

# the _g at the end of these variables helps me remember they are global
my $output_root_g;
my $XML_file_g = '';     #new global to store the loaded file   #hc
my $parallelize_bool_g = "Y";
my $fwd_primer_g = "GTGCCAGC[AC]GCCGCGGTAA";
my $rev_primer_g = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
my $seq_type_g = "PE w/ Overlap";
my $fwd_mer_len_g = 8;
my $rev_mer_len_g = 5;
my $fwd_max_shifts_g = 5;
my $rev_max_shifts_g = 5;
my $fwd_linker_g = "GA";
my $rev_linker_g = "AC";
my $min_con_depth_g = 2;
my $diginorm_max_g = "NA";
my $keep_tmp_files_g = "N";
my $con_algo_g = "NoMSA";

# flash variables
my $m_g = 30;
my $M_g = 250;
my $x_g = 0.25;
my $p_g = 33;
my $r_g = 250;
my $fr_g = 310;  # NOTE this is $fr because the frame variable is $f
my $s_g = 20;

# QC variables
my $con_trim_to_base_g = "NA";
my $con_min_len_g = 1;
my $con_min_avg_qual_g = 0;
my $con_allow_ambig_g = "N";
my $con_min_c_score_g = 0;
my $SRC_trim_to_base_g = "NA";
my $SRC_min_len_g = 1;
my $SRC_min_avg_qual_g = 0;
my $SRC_allow_gaps_g = "N";
my $SRC_allow_ambig_g = "N";


################################################################################
# Build the Main Window
my $mw = MainWindow->new();
$mw->title("MTToolbox");

# Create the notebook and fill the whole window
my $nb = $mw->NoteBook()->pack(-expand => 1,
                               -fill => 'both',
                               );


### Page 1 -- General Parameters
my $p1 = $nb->add('page1', -label => 'General');

# Page 1 fields
_output_root_field(\$output_root_g, $p1);
_parallelize_field(\$parallelize_bool_g, $p1);
_fwd_primer_field(\$fwd_primer_g, $p1);
_rev_primer_field(\$rev_primer_g, $p1);
_seq_type_field(\$seq_type_g, $p1);
_keep_tmp_files_field(\$keep_tmp_files_g, $p1);

# a padding label
{
    my $f = $p1->Frame();
    $f->Label(-text => "")->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

# Build the list box
{
    my $f = $p1->Frame;
    $list = $f->Scrolled(qw/Listbox -setgrid 1 -height 15 -scrollbars e -width 45/);
    $list->pack(-side => 'left',
                -expand => 0,
                -fill => 'y',
                -anchor => 'w');
    $list->focus;
    $list->activate(0);
    
    # Build the add button
    $f->Button(-text => 'Add Sample',
                -command => sub{add_sample($p1, $list)}
                )->pack(-side => 'top',
                        -anchor => 'nw');
    
    # Build the remove button
    $f->Button(-text => 'Remove Sample',
                -command => sub{remove_sample($list)}
                )->pack(-side => 'top',
                        -anchor => 'nw');
    
    # Build the add sample batch button
    $f->Button(-text => 'Add Sample Batch',
                -command => sub{add_sample_batch($list)}
                )->pack(-side => 'top',
                        -anchor => 'nw');
    
    # set up double click option for added samples
    $list -> bind ('<Double-ButtonPress-1>', sub { double_click($f, $list) });
    
    $f->pack(-side => 'left');
}




### Page 2 -- Advanced parameters
my $p2 = $nb->add('page2', -label => 'Advanced');

# Page 2 fields
_fwd_mer_len_field(\$fwd_mer_len_g, $p2);
_rev_mer_len_field(\$rev_mer_len_g, $p2);
_fwd_max_shifts_field(\$fwd_max_shifts_g, $p2);
_rev_max_shifts_field(\$rev_max_shifts_g, $p2);
_fwd_linker_field(\$fwd_linker_g, $p2);
_rev_linker_field(\$rev_linker_g, $p2);
_min_con_depth_field(\$min_con_depth_g, $p2);
_diginorm_max_field(\$diginorm_max_g, $p2);
_con_algo_field(\$con_algo_g, $p2);


### Page 3 -- FLASH parameters
my $p3 = $nb->add('page3', -label => 'FLASH');

_min_overlap_field(\$m_g, $p3);
_max_overlap_field(\$M_g, $p3);
_mismatch_ratio_field(\$x_g, $p3);
_phred_offset_field(\$p_g, $p3);
_avg_read_len_field(\$r_g, $p3);
_avg_frag_len_field(\$fr_g, $p3);
_sd_field(\$s_g, $p3);


### Page 4 - Quality Control Parameters
my $p4 = $nb->add('page4', -label => 'QC');

## Consensus seq QC
{
    my $f = $p4->Frame();
    my $lab = $f->Label(-text => "Conensus Seqs",
                        -font => 'arialblack 13')->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

_con_trim_to_base_field(\$con_trim_to_base_g, $p4);
_con_min_len_field(\$con_min_len_g, $p4);
_con_min_avg_qual_field(\$con_min_avg_qual_g, $p4);
_con_allow_ambig_field(\$con_allow_ambig_g, $p4);
_con_min_c_score_field(\$con_min_c_score_g, $p4);

# a padding label
{
    my $f = $p4->Frame();
    $f->Label(-text => "")->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

## Single read categories (SRCs) QC
{
    my $f = $p4->Frame();
    my $lab = $f->Label(-text => "Single Read Categories",
                        -font => 'arialblack 13')->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

_SRC_trim_to_base_field(\$SRC_trim_to_base_g, $p4);
_SRC_min_len_field(\$SRC_min_len_g, $p4);
_SRC_min_avg_qual_field(\$SRC_min_avg_qual_g, $p4);
_SRC_allow_gaps_field(\$SRC_allow_gaps_g, $p4);
_SRC_allow_ambig_field(\$SRC_allow_ambig_g, $p4);

# a padding label
{
    my $f = $p4->Frame();
    $f->Label(-text => "")->pack(-side => 'left');
    $f->pack(-fill => 'x');
}


_load_XML_GUI(\$XML_file_g, $mw);       #call to add the buttons to handle loading an XML file  #hc

# add the run button to the bottom
$mw->Button(-text => "Run",
            -command => \&run,
            )->pack(-side => 'left',
                    );

# add the exit button to the bottom
$mw->Button(-text => "Exit",
            -command => sub{exit},
            )->pack(-side => 'left',
                    );

# add a button to save current config # hc
$mw->Button(-text => "Save Current Config",
            -command => \&save_current_config,
            )->pack(-side => 'right');

MainLoop;




################################################################################
# Subroutines
sub run {
    print "running\n";
    print "Output Root: $output_root_g\n";
    print "Parallelize: $parallelize_bool_g\n";
    print "Seq Type: $seq_type_g\n";
    print "Keep Tmp Files: $keep_tmp_files_g\n";
    print "Fwd Primer: $fwd_primer_g\n";
    print "Rev Primer: $rev_primer_g\n";
    print "Fwd Linker: $fwd_linker_g\n";
    print "Rev Linker: $rev_linker_g\n";
    print "Fwd Mer len: $fwd_mer_len_g\n";
    print "Rev Mer Len: $rev_mer_len_g\n";
    print "fwd max shifts: $fwd_max_shifts_g\n";
    print "rev max shifts: $rev_max_shifts_g\n";
    print "Min Con Depth: $min_con_depth_g\n";
    print "Diginorm Max: $diginorm_max_g\n";
    print "Consensus Algorithm: $con_algo_g\n";
    
    print "\n";
    print "FLASH Params\n";
    print "m: $m_g\n";
    print "M: $M_g\n";
    print "x: $x_g\n";
    print "p: $p_g\n";
    print "r: $r_g\n";
    print "f: $fr_g\n";
    print "s: $s_g\n";
    
    print "\n";
    print "QC Paremeters\n";
    print "con_trim_to_base: $con_trim_to_base_g\n";
    print "con_min_len: $con_min_len_g\n";
    print "con_avg_qual: $con_min_avg_qual_g\n";
    print "con_allow_ambig: $con_allow_ambig_g\n";
    print "con_min_c_score: $con_min_c_score_g\n";
    print "SRC_trim_to_base: $SRC_trim_to_base_g\n";
    print "SRC_min_len: $SRC_min_len_g\n";
    print "SRC_min_avg_qual: $SRC_min_avg_qual_g\n";
    print "SRC_allow_gaps: $SRC_allow_gaps_g\n";
    print "SRC_allow_ambig: $SRC_allow_ambig_g\n";
    
    # translate type to seq_type coding
    my $seq_type_code = _translate_seq_type($seq_type_g);
    
    # set up the xml data hash ref    
    my $xml_href = {
        output_root => $output_root_g,
		split_samples_root_dir => '',
        parallelize => $parallelize_bool_g,
        fwd_primer => $fwd_primer_g,
        rev_primer => $rev_primer_g,
        seq_type => $seq_type_code,
        keep_tmp_files => $keep_tmp_files_g,
        fwd_linker => $fwd_linker_g,
        rev_linker => $rev_linker_g,
        fwd_mer_len => $fwd_mer_len_g,
        rev_mer_len => $rev_mer_len_g,
        fwd_max_shifts => $fwd_max_shifts_g,
        rev_max_shifts => $rev_max_shifts_g,
        min_con_depth => $min_con_depth_g,
        diginorm_max => $diginorm_max_g,
        con_algo => $con_algo_g,
        flash_params => {
            m => $m_g,
            M => $M_g,
            x => $x_g,
            p => $p_g,
            r => $r_g,
            f => $fr_g,
            s => $s_g,
        },
        qc_params => {
            con_trim_to_base => $con_trim_to_base_g,
            con_min_len => $con_min_len_g,
            con_min_avg_qual => $con_min_avg_qual_g,
            con_allow_ambig => $con_allow_ambig_g,
            con_min_c_score => $con_min_c_score_g,
            SRC_trim_to_base => $SRC_trim_to_base_g,
            SRC_min_len => $SRC_min_len_g,
            SRC_min_avg_qual => $SRC_min_avg_qual_g,
            SRC_allow_gaps => $SRC_allow_gaps_g,
            SRC_allow_ambig => $SRC_allow_ambig_g,
        },
        sample => \@sample_arr,
    };
    
    # Build a params object
    my $params_obj = MTToolbox::MTParams->new({href => $xml_href});
    
    # check the paramters
    my $continue_bool = 1;
    eval{
        $params_obj->check_params();
    };
    if ( $@ ) {
        $continue_bool = _handle_errors($@);
    }
    
    # If I find some error where I don't want to conintue then stop!
    if ( ! $continue_bool ) {return;} 
    
    # check the environment
    my $mt_env = MTToolbox::MTEnv->new();
    eval{
        $mt_env->check_env();
    };
    if ( $@ ) {
        $continue_bool = _handle_errors($@);
    }
    
    # If I find some error where I don't want to conintue then stop!
    if ( ! $continue_bool ) {return;}    

    # print the config file
    my $config_file = $output_root_g . "/config.xml";
    $params_obj->print_xml_file($config_file);
    
    # Run MTDriver.pl
    if ( $parallelize_bool_g =~ m/Y/i ) {
        my $command = "bsub -q week -n 1 " .
                      "-o $output_root_g/MTDriver.log " .
                      "-e $output_root_g/MTDriver.err " .
                      "-J MTDriver " .
                      "MTDriver.pl --params_file $config_file";
        system("$command");
        
        # Tell the user it is running and ask if they want to run another
        my $ans = $mw->messageBox(-title => "Status",
                                  -message => "MTToolbox is running.\n\n" .
                                              "Start another run?",
                                  -type => "YesNo",
                                  -default => 'No');
        if ( $ans =~ m/Yes/i ) {
            reset_GUI();
        }
        else {
            exit;
        }
    }
    else {
        system("MTDriver.pl --params_file $config_file");
    }
}

sub add_sample {
    my $w = shift;
    my $list = shift;
    
    # Variables for a sample
    my $id;
    my $fwd_file;
    my $rev_file;
    my $barcode;
    my $seq_type_l = $seq_type_g;
    my $fwd_primer_l = $fwd_primer_g;
    my $rev_primer_l = $rev_primer_g;
    my $fwd_mer_len_l = $fwd_mer_len_g;
    my $rev_mer_len_l = $rev_mer_len_g;
    my $fwd_max_shifts_l = $fwd_max_shifts_g;
    my $rev_max_shifts_l = $rev_max_shifts_g;
    my $fwd_linker_l = $fwd_linker_g;
    my $rev_linker_l = $rev_linker_g;
    my $min_con_depth_l = $min_con_depth_g;
    my $diginorm_max_l = $diginorm_max_g;
    my $con_algo_l = $con_algo_g;
    my $m_l = $m_g;
    my $M_l = $M_g;
    my $x_l = $x_g;
    my $p_l = $p_g;
    my $r_l = $r_g;
    my $fr_l = $fr_g;  # NOTE this is $fr because the frame variable is $f
    my $s_l = $s_g;
    my $con_trim_to_base_l = $con_trim_to_base_g;
    my $con_min_len_l = $con_min_len_g;
    my $con_min_avg_qual_l = $con_min_avg_qual_g;
    my $con_allow_ambig_l = $con_allow_ambig_g;
    my $con_min_c_score_l = $con_min_c_score_g;
    my $SRC_trim_to_base_l = $SRC_trim_to_base_g;
    my $SRC_min_len_l = $SRC_min_len_g;
    my $SRC_min_avg_qual_l = $SRC_min_avg_qual_g;
    my $SRC_allow_gaps_l = $SRC_allow_gaps_g;
    my $SRC_allow_ambig_l = $SRC_allow_ambig_g;
    
    # Create the dialog box
    $id = "P" . $list->size();
    my $db = $w->DialogBox(-title => $id,
                           -buttons => ["Add", "Cancel"],
                           -default_button => "Add");
    
    # everything in the dialog box is in a notebook object
    my $nb = $db->NoteBook()->pack(-expand => 1,
                                   -fill => 'both',
                                );
    
    # Page 1 -- General Parameters
    my $p1 = $nb->add('page1', -label => 'General');
    
    _fwd_read_file_field(\$fwd_file, $p1);
    _rev_read_file_field(\$rev_file, $p1);
    _seq_type_field(\$seq_type_l, $p1);
    _fwd_primer_field(\$fwd_primer_l, $p1);
    _rev_primer_field(\$rev_primer_l, $p1);

    
    # Page 2 -- Advanced Parameters
    my $p2 = $nb->add('page2', -label => 'Advanced');
    
    _fwd_mer_len_field(\$fwd_mer_len_l, $p2); 
    _rev_mer_len_field(\$rev_mer_len_l, $p2); 
    _fwd_max_shifts_field(\$fwd_max_shifts_l, $p2);
    _rev_max_shifts_field(\$rev_max_shifts_l, $p2);
    _fwd_linker_field(\$fwd_linker_l, $p2);
    _rev_linker_field(\$rev_linker_l, $p2);
    _min_con_depth_field(\$min_con_depth_l, $p2);
    _diginorm_max_field(\$diginorm_max_l, $p2);
    _con_algo_field(\$con_algo_l, $p2);
    
    
    # Page 3 -- FLASH parameters
    my $p3 = $nb->add('page3', -label => 'FLASH');
    
    _min_overlap_field(\$m_l, $p3);
    _max_overlap_field(\$M_l, $p3);
    _mismatch_ratio_field(\$x_l, $p3);
    _phred_offset_field(\$p_l, $p3);
    _avg_read_len_field(\$r_l, $p3);
    _avg_frag_len_field(\$fr_l, $p3);
    _sd_field(\$s_l, $p3);
    
    # Page 4 -- QC parameters
    my $p4 = $nb->add('page4', -label => 'QC');

    # Consensus seq QC
    {
        my $f = $p4->Frame();
        my $lab = $f->Label(-text => "Conensus Seqs",
                            -font => 'arialblack 13')->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    _con_trim_to_base_field(\$con_trim_to_base_l, $p4);
    _con_min_len_field(\$con_min_len_l, $p4);
    _con_min_avg_qual_field(\$con_min_avg_qual_l, $p4);
    _con_allow_ambig_field(\$con_allow_ambig_l, $p4);
    _con_min_c_score_field(\$con_min_c_score_l, $p4);
    
    # a padding label
    {
        my $f = $p4->Frame();
        $f->Label(-text => "")->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    
    # Single read categories (SRCs) QC
    {
        my $f = $p4->Frame();
        my $lab = $f->Label(-text => "Single Read Categories",
                            -font => 'arialblack 13')->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    _SRC_trim_to_base_field(\$SRC_trim_to_base_l, $p4);
    _SRC_min_len_field(\$SRC_min_len_l, $p4);
    _SRC_min_avg_qual_field(\$SRC_min_avg_qual_l, $p4);
    _SRC_allow_gaps_field(\$SRC_allow_gaps_l, $p4);
    _SRC_allow_ambig_field(\$SRC_allow_ambig_l, $p4);

    # a padding label
    {
        my $f = $p4->Frame();
        $f->Label(-text => "")->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    
    
    # listen for add button
    my $answer = $db->Show();
    if ( $answer eq 'Add' ) {
        # Guess the name based on the file
        my $name = guess_file_name($fwd_file);
        print "name: $name\n";
        
        # add to listbox
        $list->insert('end', "$name");
        $list->pack();
        
        # Make the sample hash
        my $sample_hash = {sample_id => $id,
                            fwd_file => $fwd_file,
                            rev_file => $rev_file,
                            barcode => $barcode};
        
        # The local parameters I only add if they are differen than the
        # corresponding global parameter (i.e. it has been updated)
        if ( $seq_type_l ne $seq_type_g ) {
            $sample_hash->{seq_type} = $seq_type_l;
        }
        if ( $fwd_primer_l ne $fwd_primer_g ) {
            $sample_hash->{fwd_primer} = $fwd_primer_l;
        }
        if ( $rev_primer_l ne $rev_primer_g ) {
            $sample_hash->{rev_primer} = $rev_primer_l;
        }
        if ( $fwd_mer_len_l ne $fwd_mer_len_g ) {
            $sample_hash->{fwd_mer_len} = $fwd_mer_len_l;
        }
        if ( $rev_mer_len_l ne $rev_mer_len_g ) {
            $sample_hash->{rev_mer_len} = $rev_mer_len_l;
        }
        if ( $fwd_max_shifts_l ne $fwd_max_shifts_g ) {
            $sample_hash->{fwd_max_shifts} = $fwd_max_shifts_l;
        }
        if ( $rev_max_shifts_l ne $rev_max_shifts_g ) {
            $sample_hash->{rev_max_shifts} = $rev_max_shifts_l;
        }
        if( $fwd_linker_l ne $fwd_linker_g ) {
            $sample_hash->{fwd_linker} = $fwd_linker_l;
        }
        if ( $rev_linker_l ne $rev_linker_g ) {
            $sample_hash->{rev_linker} = $rev_linker_l;
        }
        if ( $min_con_depth_l ne $min_con_depth_g ) {
            $sample_hash->{min_con_depth} = $min_con_depth_l;
        }
        if ( $diginorm_max_l ne $diginorm_max_g ) {
            $sample_hash->{diginorm_max} = $diginorm_max_l;
        }
        if ( $con_algo_l ne $con_algo_g ) {
            $sample_hash->{con_algo} = $con_algo_l;
        }
        if ( $m_l ne $m_g ) {
            $sample_hash->{flash_params}->{m} = $m_l;
        }
        if ( $M_l ne $M_g ) {
            $sample_hash->{flash_params}->{M} = $M_l;
        }
        if ( $x_l ne $x_g ) {
            $sample_hash->{flash_params}->{x} = $x_l;
        }
        if ( $p_l ne $p_g ) {
            $sample_hash->{flash_params}->{p} = $p_l;
        }
        if ( $r_l ne $r_g ) {
            $sample_hash->{flash_params}->{r} = $r_l;
        }
        if ( $fr_l ne $fr_g ) {
            $sample_hash->{flash_params}->{f} = $fr_l;
        }
        if ( $s_l ne $s_g ) {
            $sample_hash->{flash_params}->{s} = $s_l;
        }
        if ( $con_trim_to_base_l ne $con_trim_to_base_g ) {
            $sample_hash->{qc_params}->{con_trim_to_base} = $con_trim_to_base_l;
        }
        if ( $con_min_len_l ne $con_min_len_g ) {
            $sample_hash->{qc_params}->{con_min_len} = $con_min_len_l;
        }
        if ( $con_min_avg_qual_l ne $con_min_avg_qual_g ) {
            $sample_hash->{qc_params}->{con_min_avg_qual} = $con_min_avg_qual_l;
        }
        if ( $con_allow_ambig_l ne $con_allow_ambig_g ) {
            $sample_hash->{qc_params}->{con_allow_ambig} = $con_allow_ambig_l;
        }
        if ( $con_min_c_score_l ne $con_min_c_score_g ) {
            $sample_hash->{qc_params}->{con_min_c_score} = $con_min_c_score_l;
        }
        if ( $SRC_trim_to_base_l ne $SRC_trim_to_base_g ) {
            $sample_hash->{qc_params}->{SRC_trim_to_base} = $SRC_trim_to_base_l;
        }
        if ( $SRC_min_len_l ne $SRC_min_len_g ) {
            $sample_hash->{qc_params}->{SRC_min_len} = $SRC_min_len_l;
        }
        if ( $SRC_min_avg_qual_l ne $SRC_min_avg_qual_g ) {
            $sample_hash->{qc_params}->{SRC_min_avg_qual} = $SRC_min_avg_qual_l;
        }
        if ( $SRC_allow_gaps_l ne $SRC_allow_gaps_g ) {
            $sample_hash->{qc_params}->{SRC_allow_gaps} = $SRC_allow_gaps_l;
        }
        if ( $SRC_allow_ambig_l ne $SRC_allow_ambig_g ) {
            $sample_hash->{qc_params}->{SRC_allow_ambig} = $SRC_allow_ambig_l;
        }
        
        # add to sample href
        # NOTE: I need to make sure the id is unique!!!
        push @sample_arr, $sample_hash;
    }
    
    #print "answer: $answer\n";
}

sub remove_sample {
    my $list = shift;
    
    # remove from the listbox
    my @selected_i = sort $list->curselection();
    $list->delete($selected_i[0], $selected_i[-1]);
    $list->pack();
    
    # remove from the sample data structure
    foreach my $i ( @selected_i ) {
        splice @sample_arr, $i, 1;
    }
    
    recalculate_p_nums();
}

sub add_sample_batch {
    my ($list) = @_;
    
    my $dir;  # a directory to look in for the sample batch
    
    # NOTE: I have to make this frame object in order to use dirDialog.  However,
    # I don't actually add it to the GUI (e.g. I don't call pack).
    my $f = $mw->Frame();
    my $ent = $f->Entry(-width => 20,
                        -textvariable => \$dir);
    
    dirDialog($f, $ent);
    
    opendir(DIR, $dir) or croak("Cannot open dir: $dir");
    
    my %fwd_files;
    my %rev_files;
    my %weird_files;
    FILE: while (my $file = readdir(DIR)) {
        print "file: $file\n";
        
        # skip files that begin with a dot
        if ( $file =~ m/^\./ ) {next FILE;}

        if ( $file =~ m/.*\_R1\_.*/ ) {
            $fwd_files{$file} = 1;
        }
        if ( $file =~ m/.*_R2_.*/ ) {
            $rev_files{$file} = 1;
        }
        else {
            $weird_files{$file} = 1;
        }
    }
        
    # Find each fwd file match in the rev_file hash
    my $sample_count = $list->size();
    foreach my $fwd_file ( keys %fwd_files ) {
        my $rev_file_temp = $fwd_file;
        $rev_file_temp =~ s/_R1_/_R2_/;
        
        # add this sample to the list box ($list)
        my $name = guess_file_name($fwd_file);
        $list->insert('end', "$name");
        $list->pack();
        
        # add this sample to the sample list (@sample_arr)
        my $id = "P" . $sample_count++;
        print "id: $id\n";
        
        # if there is a rev file match then it must be a PE sample
        if ( defined $rev_files{$rev_file_temp} ) {
            push @sample_arr, {sample_id => $id,
                               fwd_file => "$dir/$fwd_file",
                               rev_file => "$dir/$rev_file_temp",
                               barcode => guess_barcode($fwd_file)};
            
            # delete that rev file
            delete $rev_files{$rev_file_temp};
        }
        else {
            push @sample_arr, {sample_id => $id,
                               fwd_file => "$dir/$fwd_file",
                               barcode => guess_barcode($fwd_file)};
        }
        
        # delete that fwd file.  rev files are deleted in the if statement
        delete $fwd_files{$fwd_file};
    }
    
    close(DIR);
}

sub double_click {
    my ($w, $list) = @_;
    
    my $selected_i = @{$list->curselection()}[0];
    print "Selected: $selected_i\n";
    my $sample_hash = $sample_arr[$selected_i];
    
    # Variables for this sample
    my $id = $sample_hash->{sample_id};
    my $fwd_file = $sample_hash->{fwd_file};
    my $rev_file = $sample_hash->{rev_file};
    my $barcode = $sample_hash->{barcode};
    my $seq_type_l = $sample_hash->{seq_type};
    my $fwd_primer_l = $sample_hash->{fwd_primer};
    my $rev_primer_l = $sample_hash->{rev_primer};
    my $fwd_mer_len_l = $sample_hash->{fwd_mer_len};
    my $rev_mer_len_l = $sample_hash->{rev_mer_len};
    my $fwd_max_shifts_l = $sample_hash->{fwd_max_shifts};
    my $rev_max_shifts_l = $sample_hash->{rev_max_shifts};
    my $fwd_linker_l = $sample_hash->{fwd_linker};
    my $rev_linker_l = $sample_hash->{rev_linker};
    my $min_con_depth_l = $sample_hash->{min_con_depth};
    my $diginorm_max_l = $sample_hash->{diginorm_max};
    my $con_algo_l = $sample_hash->{con_algo};
    my $m_l = $sample_hash->{flash_params}->{m};
    my $M_l = $sample_hash->{flash_params}->{M};
    my $x_l = $sample_hash->{flash_params}->{x};
    my $p_l = $sample_hash->{flash_params}->{p};
    my $r_l = $sample_hash->{flash_params}->{r};
    my $fr_l = $sample_hash->{flash_params}->{f};
    my $s_l = $sample_hash->{flash_params}->{s};
    my $con_trim_to_base_l = $sample_hash->{qc_params}->{con_trim_to_base};
    my $con_min_len_l = $sample_hash->{qc_params}->{con_min_len};
    my $con_min_avg_qual_l = $sample_hash->{qc_params}->{con_min_avg_qual};
    my $con_allow_ambig_l = $sample_hash->{qc_params}->{con_allow_ambig};
    my $con_min_c_score_l = $sample_hash->{qc_params}->{con_min_c_score};
    my $SRC_trim_to_base_l = $sample_hash->{qc_params}->{SRC_trim_to_base};
    my $SRC_min_len_l = $sample_hash->{qc_params}->{SRC_min_len};
    my $SRC_min_avg_qual_l = $sample_hash->{qc_params}->{SRC_min_avg_qual};
    my $SRC_allow_gaps_l = $sample_hash->{qc_params}->{SRC_allow_gaps};
    my $SRC_allow_ambig_l = $sample_hash->{qc_params}->{SRC_allow_ambig};

    
    # Any variable that are not defined locally in the sample default to global
    # variables.
    
    # General variables
    if ( ! defined $seq_type_l ) { $seq_type_l = $seq_type_g; }
    if ( ! defined $fwd_primer_l ) { $fwd_primer_l = $fwd_primer_g; }
    if ( ! defined $rev_primer_l ) { $rev_primer_l = $rev_primer_g; }
    
    # Advanced variables
    if ( ! defined $fwd_mer_len_l ) { $fwd_mer_len_l = $fwd_mer_len_g; }
    if ( ! defined $rev_mer_len_l ) { $rev_mer_len_l = $rev_mer_len_g; }
    if ( ! defined $fwd_max_shifts_l ) { $fwd_max_shifts_l = $fwd_max_shifts_g; }
    if ( ! defined $rev_max_shifts_l ) { $rev_max_shifts_l = $rev_max_shifts_g; }
    if ( ! defined $fwd_linker_l ) { $fwd_linker_l = $fwd_linker_g; }
    if ( ! defined $rev_linker_l ) { $rev_linker_l = $rev_linker_g; }
    if ( ! defined $min_con_depth_l ) { $min_con_depth_l = $min_con_depth_g; }
    if ( ! defined $diginorm_max_l ) { $diginorm_max_l = $diginorm_max_g; }
    if ( ! defined $con_algo_l ) { $con_algo_l = $con_algo_g; }
    
    # flash variables
    if ( ! defined $m_l ) { $m_l = $m_g; }
    if ( ! defined $M_l ) { $M_l = $M_g; }
    if ( ! defined $x_l ) { $x_l = $x_g; }
    if ( ! defined $p_l ) { $p_l = $p_g; }
    if ( ! defined $r_l ) { $r_l = $r_g; }
    if ( ! defined $fr_l ) { $fr_l = $fr_g; }
    if ( ! defined $s_l ) { $s_l = $s_g; }
    
    # Quality control variables
    if ( ! defined $con_trim_to_base_l ) { $con_trim_to_base_l = $con_trim_to_base_g; }
    if ( ! defined $con_min_len_l ) { $con_min_len_l = $con_min_len_g; }
    if ( ! defined $con_min_avg_qual_l ) { $con_min_avg_qual_l = $con_min_avg_qual_g; }
    if ( ! defined $con_allow_ambig_l ) { $con_allow_ambig_l = $con_allow_ambig_g; }
    if ( ! defined $con_min_c_score_l ) { $con_min_c_score_l = $con_min_c_score_g; }
    if ( ! defined $SRC_trim_to_base_l ) { $SRC_trim_to_base_l = $SRC_trim_to_base_g; }
    if ( ! defined $SRC_min_len_l ) { $SRC_min_len_l = $SRC_min_len_g; }
    if ( ! defined $SRC_min_avg_qual_l ) { $SRC_min_avg_qual_l = $SRC_min_avg_qual_g; }
    if ( ! defined $SRC_allow_gaps_l ) { $SRC_allow_gaps_l = $SRC_allow_gaps_g; }
    if ( ! defined $SRC_allow_ambig_l ) { $SRC_allow_ambig_l = $SRC_allow_ambig_g; }

    
    
    # add the dialog box
    my $db = $w->DialogBox(-title => $id,
                           -buttons => ["Update", "Cancel"],
                           -default_button => "Update");
    
    # everything in the dialog box is in a notebook object
    my $nb = $db->NoteBook()->pack(-expand => 1,
                                   -fill => 'both',
                                );
    
    # Page 1 -- General Parameters
    my $p1 = $nb->add('page1', -label => 'General');
    
    _fwd_read_file_field(\$fwd_file, $p1);
    _rev_read_file_field(\$rev_file, $p1);
    _seq_type_field(\$seq_type_l, $p1);
    _fwd_primer_field(\$fwd_primer_l, $p1);
    _rev_primer_field(\$rev_primer_l, $p1);
    
    # Page 2 -- Advanced parameters
    my $p2 = $nb->add('page2', -label => 'Advanced');
    
    _fwd_mer_len_field(\$fwd_mer_len_l, $p2); 
    _rev_mer_len_field(\$rev_mer_len_l, $p2); 
    _fwd_max_shifts_field(\$fwd_max_shifts_l, $p2);
    _rev_max_shifts_field(\$rev_max_shifts_l, $p2);
    _fwd_linker_field(\$fwd_linker_l, $p2);
    _rev_linker_field(\$rev_linker_l, $p2);
    _min_con_depth_field(\$min_con_depth_l, $p2);
    _diginorm_max_field(\$diginorm_max_l, $p2);
    _con_algo_field(\$con_algo_l, $p2);
    
    
    # Page 3 -- FLASH parameters
    my $p3 = $nb->add('page3', -label => 'FLASH');
    
    _min_overlap_field(\$m_l, $p3);
    _max_overlap_field(\$M_l, $p3);
    _mismatch_ratio_field(\$x_l, $p3);
    _phred_offset_field(\$p_l, $p3);
    _avg_read_len_field(\$r_l, $p3);
    _avg_frag_len_field(\$fr_l, $p3);
    _sd_field(\$s_l, $p3);
    
    # Page 4 -- Quality Control Parameters
    my $p4 = $nb->add('page4', -label => 'QC');

    # Consensus seq QC
    {
        my $f = $p4->Frame();
        my $lab = $f->Label(-text => "Conensus Seqs",
                            -font => 'arialblack 13')->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    _con_trim_to_base_field(\$con_trim_to_base_l, $p4);
    _con_min_len_field(\$con_min_len_l, $p4);
    _con_min_avg_qual_field(\$con_min_avg_qual_l, $p4);
    _con_allow_ambig_field(\$con_allow_ambig_l, $p4);
    _con_min_c_score_field(\$con_min_c_score_l, $p4);
    
    # a padding label
    {
        my $f = $p4->Frame();
        $f->Label(-text => "")->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    
    # Single read categories (SRCs) QC
    {
        my $f = $p4->Frame();
        my $lab = $f->Label(-text => "Single Read Categories",
                            -font => 'arialblack 13')->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }
    _SRC_trim_to_base_field(\$SRC_trim_to_base_l, $p4);
    _SRC_min_len_field(\$SRC_min_len_l, $p4);
    _SRC_min_avg_qual_field(\$SRC_min_avg_qual_l, $p4);
    _SRC_allow_gaps_field(\$SRC_allow_gaps_l, $p4);
    _SRC_allow_ambig_field(\$SRC_allow_ambig_l, $p4);
    
    # a padding label
    {
        my $f = $p4->Frame();
        $f->Label(-text => "")->pack(-side => 'left');
        $f->pack(-fill => 'x');
    }

    # wait for a button to be pressed
    my $answer = $db->Show();
    
    # Update the values from the user
    if ( $answer eq 'Update' ) {
        # Guess the name based on the file
        #my $name = guess_file_name($fwd_file);
        #print "name: $name\n";
        
        # Update the sample hash with everything that changes
        my $new_sample_hash = {sample_id => $id,
                               fwd_file => $fwd_file,
                               rev_file => $rev_file,
                               barcode => $barcode};
        
        if ( $seq_type_l ne $seq_type_g ) {
            $new_sample_hash->{seq_type} = $seq_type_l;
        }
        if ( $fwd_primer_l ne $fwd_primer_g ) {
            $new_sample_hash->{fwd_primer} = $fwd_primer_l;
        }
        if ( $rev_primer_l ne $rev_primer_g ) {
            $new_sample_hash->{rev_primer} = $rev_primer_l;
        }
        if ( $fwd_mer_len_l ne $rev_mer_len_g ) {
            $new_sample_hash->{fwd_mer_len} = $fwd_mer_len_l;
        }
        if ( $rev_mer_len_l ne $rev_mer_len_g ) {
            $new_sample_hash->{rev_mer_len} = $rev_mer_len_l;
        }
        if ( $fwd_max_shifts_l ne $fwd_max_shifts_g ) {
            $new_sample_hash->{fwd_max_shifts} = $fwd_max_shifts_l;
        }
        if ( $rev_max_shifts_l ne $rev_max_shifts_g ) {
            $new_sample_hash->{rev_max_shifts} = $rev_max_shifts_l;
        }
        if( $fwd_linker_l ne $fwd_linker_g ) {
            $new_sample_hash->{fwd_linker} = $fwd_linker_l;
        }
        if ( $rev_linker_l ne $rev_linker_g ) {
            $new_sample_hash->{rev_linker} = $rev_linker_l;
        }
        if ( $min_con_depth_l ne $min_con_depth_g ) {
            $new_sample_hash->{min_con_depth} = $min_con_depth_l;
        }
        if ( $diginorm_max_l ne $diginorm_max_g ) {
            $new_sample_hash->{diginorm_max} = $diginorm_max_l;
        }
        if ( $con_algo_l ne $con_algo_g ) {
            $new_sample_hash->{con_algo} = $con_algo_l;
        }
        
        # update the flash parameters
        my $flash_params_href = {};
        if ( $m_l ne $m_g ) {
            $flash_params_href->{m} = $m_l;
        }
        if ( $M_l ne $M_g ) {
            $flash_params_href->{M} = $M_l;
        }
        if ( $x_l ne $x_g ) {
            $flash_params_href->{x} = $x_l;
        }
        if ( $p_l ne $p_g ) {
            $flash_params_href->{p} = $p_l;
        }
        if ( $r_l ne $r_g ) {
            $flash_params_href->{r} = $r_l;
        }
        if ( $fr_l ne $fr_g ) {
            $flash_params_href->{f} = $fr_l;
        }
        if ( $s_l ne $s_g ) {
            $flash_params_href->{s} = $s_l;
        }
        # add the updated flash parameters to the new_sample_hash
        $new_sample_hash->{flash_params} = $flash_params_href;
        
        
        # update the QC parameters
        my $qc_params_href = {};
        if ( $con_trim_to_base_l ne $con_trim_to_base_g ) {
            $qc_params_href->{con_trim_to_base} = $con_trim_to_base_l;
        }
        if ( $con_min_len_l ne $con_min_len_g ) {
            $qc_params_href->{con_min_len} = $con_min_len_l;
        }
        if ( $con_min_avg_qual_l ne $con_min_avg_qual_g ) {
            $qc_params_href->{con_min_avg_qual} = $con_min_avg_qual_l;
        }
        if ( $con_allow_ambig_l ne $con_allow_ambig_g ) {
            $qc_params_href->{con_allow_ambig} = $con_allow_ambig_l;
        }
        if ( $con_min_c_score_l ne $con_min_c_score_g ) {
            $qc_params_href->{con_min_c_score} = $con_min_c_score_l;
        }
        if ( $SRC_trim_to_base_l ne $SRC_trim_to_base_g ) {
            $qc_params_href->{SRC_trim_to_base} = $SRC_trim_to_base_l;
        }
        if ( $SRC_min_len_l ne $SRC_min_len_g ) {
            $qc_params_href->{SRC_min_len} = $SRC_min_len_l;
        }
        if ( $SRC_min_avg_qual_l ne $SRC_min_avg_qual_g ) {
            $qc_params_href->{SRC_min_avg_qual} = $SRC_min_avg_qual_l;
        }
        if ( $SRC_allow_gaps_l ne $SRC_allow_gaps_g ) {
            $qc_params_href->{SRC_allow_gaps} = $SRC_allow_gaps_l;
        }
        if ( $SRC_allow_ambig_l ne $SRC_allow_ambig_g ) {
            $qc_params_href->{SRC_allow_ambig} = $SRC_allow_ambig_l;
        }
        # add the updated qc parameters ot the new_sample_hash
        $new_sample_hash->{qc_params} = $qc_params_href;
        
        
        # check for errors        
        my $continue_bool = 1;
        eval{
            MTToolbox::MTParams::check_sample_entry($new_sample_hash,
                                          $seq_type_g,
                                          '');
        };
        if ( $@ ) {
            $continue_bool = _handle_errors($@);
        }
        
        # If I find some error where I don't want to conintue then stop!
        # The nothing from the pending update will be saved if we return at
        # this point.
        if ( ! $continue_bool ) {
            return;
        }  
        
        # update the list object with the new_sample_hash
        $sample_arr[$selected_i] = $new_sample_hash;
    }
}

sub reset_GUI {
    $list->delete(0, $list->size()); # a list object that goes on the GUI.
    @sample_arr = ();  # stores all info for each sample in array
    
    # the _g at the end of these variables helps me remember they are global
    $output_root_g = '';
    $parallelize_bool_g = "Y";
    $fwd_primer_g = "GTGCCAGC[AC]GCCGCGGTAA";
    $rev_primer_g = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
    $seq_type_g = "PE w/ Overlap";
    $keep_tmp_files_g = "N";
    $fwd_mer_len_g = 8;
    $rev_mer_len_g = 5;
    $fwd_max_shifts_g = 5;
    $rev_max_shifts_g = 5;
    $fwd_linker_g = "GA";
    $rev_linker_g = "AC";
    $min_con_depth_g = 2;
    $diginorm_max_g = "NA";
    $con_algo_g = "NoMSA";
    
    # flash variables
    $m_g = 30;
    $M_g = 250;
    $x_g = 0.25;
    $p_g = 33;
    $r_g = 250;
    $fr_g = 310;  # NOTE this is $fr because the frame variable is $f
    $s_g = 20;
    
    # QC variables
    $con_trim_to_base_g = "NA";
    $con_min_len_g = 1;
    $con_min_avg_qual_g = 0;
    $con_allow_ambig_g = "N";
    $con_min_c_score_g = 0;
    $SRC_trim_to_base_g = "NA";
    $SRC_min_len_g = 1;
    $SRC_min_avg_qual_g = 0;
    $SRC_allow_gaps_g = "N";
    $SRC_allow_ambig_g = "N";
    
    $mw->update();
}

sub guess_file_name {
    my ($file) = @_;
    
    my $name = "";
    if ( $file =~ m/(.*)\/(.*)_L001_R1_001.fastq/i ) {
        $name = $2;
    }
    elsif ( $file =~ m/(.*)_L001_R1_001.fastq/i ) {
        $name = $1;
    }
    else {
        $name = $file;
    }
    
    return $name;
}

sub guess_barcode {
    my ($file) = @_;
    
    my $barcode = "";
    if ( $file =~ m/.*_([ATCG]*)_.*/i ) {
        $barcode = $1;
    }
    
    print "Barcode Guess: $barcode\n";
    
    return $barcode;
}

sub recalculate_p_nums {
    my $index = 0;
    foreach my $sample_hash ( @sample_arr ) {
        $sample_hash->{sample_id} = "P" . $index++;
    }
    
    return 1;
}

sub _translate_seq_type {
    my ($seq_type) = @_;
    
    my $seq_type_code = $seq_type;
    if ( $seq_type =~ m/SE/ ){ $seq_type_code = 1;}
    elsif ( $seq_type =~ m/PE w\/ Overlap/ ) {$seq_type_code = 2;}
    elsif ($seq_type =~ m/PE w\/o Overlap/ ) { $seq_type_code = 3;}
    
    return $seq_type_code;
}

sub _handle_errors {
    my $continue = 1;
    
    if ( my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingRequiredTag') ) {
        my $message = $e->error . "\n" .
                      "Tag: " . $e->tag_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingTagValue') ) {
        my $message = $e->error . "\n" .
                      "Tag: " . $e->tag_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::UnrecignizedTag') ) {
        my $message = $e->error . "\n" .
                      "Tag: " . $e->tag_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingExec') ) {
        my $message = $e->error . "\n" .
                      "Program: " . $e->exec_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::LongLinker') ) {
        my $message = $e->error . "\n" .
                      "Tag: " . $e->tag_name . "\n";
        my $ans = $mw->messageBox(-title => "ERROR",
                                    -message => $message,
                                    -type => "YesNo",
                                    -icon => 'warning',
                                    -default => 'Yes');
        if ( $ans =~ m/No/i ) {
            $continue = 0;
        }
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::UnrecignizedChar') ) {
        my $message = $e->error . "\n" .
                      "Orientation: " . $e->orientation() . "\n" .
                      "Length: " . $e->length . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::BadSeqType') ) {
        my $message = $e->error . "\n" .
                      "Seq_type: " . $e->value . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
    }
    elsif ( $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingExec') ) {
        my $message = $e->error . "\n" .
                      "Program: " . $e->exec_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::DoesNotExist::File') ) {
        my $message = $e->error . "\n" .
                      "File: " . $e->file_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::DoesNotExist::Dir') ) {
        my $message = $e->error . "\n" .
                      "Dir: " . $e->dir_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit') ) {
        my $message = $e->error . "\n" .
                      "Value: " . $e->value . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall') ) {
        my $message = $e->error . "\n" .
                      "Value: " . $e->value . "\n" .
                      "MIN: " . $e->MIN . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::Digit::TooBig') ) {
        my $message = $e->error . "\n" .
                      "Value: " . $e->value . "\n" .
                      "MAX: " . $e->MAX . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::Digit::OOB') ) {
        my $message = $e->error . "\n" .
                      "Value: " . $e->value . "\n" .
                      "MIN: " . $e->MIN . "\n" .
                      "MAX: " . $e->MAX . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::File::Empty') ) {
        my $message = $e->error . "\n" .
                      "File: " . $e->file_name . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::File::BadExtension') ) {
        my $message = $e->error . "\n" .
                      "File: " . $e->file_name . "\n" .
                      "Extension: " . $e->ext . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( $e = Exception::Class->caught('MyX::Generic::BadValue') ) {
        my $message = $e->error . "\n" .
                      "Value: " . $e->value() . "\n";
        $mw->messageBox(-title => "ERROR",
                        -message => $message,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');
        $continue = 0;
    }
    elsif ( defined $e ) {
        $e->rethrow();
        $continue = 0;
    }
    
    return $continue;
}

sub dirDialog {
    my $w = shift;
    my $ent = shift;
    my $dir;
    $dir = $w->chooseDirectory(initialdir => "~");
    if ( defined $dir and $dir ne '' ) {
        # populate the Output Directory Root field in the form
        $ent->delete(0, 'end');
        $ent->insert(0, $dir);
        $ent->xview('end');
    }
}

sub fileDialog {
    my $w = shift;
    my $ent = shift;
    my $operation = shift;
    my $types;
    my $file;
    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    my @types =
      (["All files",		'*'],
	   ["Text files",           [qw/.txt .doc/]],
       ["Text files",           '',             'TEXT'],
       ["Perl Scripts",         '.pl',		'TEXT'],
       ["C Source Files",	['.c', '.h']],
       ["All Source Files",     [qw/.tcl .c .h/]],
       ["Image Files",		'.gif'],
       ["Image Files",		['.jpeg', '.jpg']],
       ["Image Files",   	'',		[qw/GIFF JPEG/]],
       ["All files",		'*']
      );
    if ($operation eq 'open') {
	$file = $w->getOpenFile(-filetypes => \@types);
    } else {
	$file = $w->getSaveFile(-filetypes => \@types,
				-initialfile => 'Untitled',
				-defaultextension => '.txt');
    }
    if (defined $file and $file ne '') {
        $ent->delete(0, 'end');
        $ent->insert(0, $file);
        $ent->xview('end');
    }
}

# <hc>
sub load_XML {  
    my $xml_hash;
    my $err_msg;

    reset_GUI();

    #catch errors while reading in the file
    eval{ $xml_hash = XMLin($XML_file_g); };
    if ( $@ ) {
        
        if ( $XML_file_g eq '' ) { $err_msg = "No XML file specified."; }
        else                     { $err_msg = "Error reading XML file: $XML_file_g"; }
        
        $mw->messageBox(-title => "ERROR",
                        -message => $err_msg,
                        -type => "OK",
                        -icon => 'error',
                        -default => 'OK');

        return
    }
    
    #hash to convert $seq_type integer in xml to str for GUI
    my %conv_seq_type = (   '1' =>  'SE',
                            '2' =>  'PE w/ Overlap',
                            '3' =>  'PE w/o Overlap'
                        );

    #set global parameters (denoted with _g)
    $output_root_g = $xml_hash->{'output_root'};
    $parallelize_bool_g = $xml_hash->{'parallelize'};
    $fwd_primer_g = $xml_hash->{'fwd_primer'};
    $rev_primer_g = $xml_hash->{'rev_primer'};
    $seq_type_g = $conv_seq_type{$xml_hash->{'seq_type'}};      #this uses the xml (integer) field as a key in the conversion hash to get the string value
    $keep_tmp_files_g = $xml_hash->{'keep_tmp_files'};
    $fwd_mer_len_g = $xml_hash->{'fwd_mer_len'};
    $rev_mer_len_g = $xml_hash->{'rev_mer_len'};
    $fwd_max_shifts_g = $xml_hash->{'fwd_max_shifts'};
    $rev_max_shifts_g = $xml_hash->{'rev_max_shifts'};
    $fwd_linker_g = $xml_hash->{'fwd_linker'};
    $rev_linker_g = $xml_hash->{'rev_linker'};
    $min_con_depth_g = $xml_hash->{'min_con_depth'};
    $diginorm_max_g = $xml_hash->{'diginorm_max'};
    $con_algo_g = $xml_hash->{'con_algo'};
    
    # flash variables
    $m_g = $xml_hash->{'flash_params'}{'m'};
    $M_g = $xml_hash->{'flash_params'}{'M'};
    $x_g = $xml_hash->{'flash_params'}{'x'};
    $p_g = $xml_hash->{'flash_params'}{'p'};
    $r_g = $xml_hash->{'flash_params'}{'r'};
    $fr_g = $xml_hash->{'flash_params'}{'f'};
    $s_g = $xml_hash->{'flash_params'}{'s'};
    
    # QC variables
    $con_trim_to_base_g = $xml_hash->{'qc_params'}{'con_trim_to_base'};
    $con_min_len_g = $xml_hash->{'qc_params'}{'con_min_len'};
    $con_min_avg_qual_g = $xml_hash->{'qc_params'}{'con_min_avg_qual'};
    $con_allow_ambig_g = $xml_hash->{'qc_params'}{'con_allow_ambig'};
    $con_min_c_score_g = $xml_hash->{'qc_params'}{'con_min_c_score'};
    $SRC_trim_to_base_g = $xml_hash->{'qc_params'}{'SRC_trim_to_base'};
    $SRC_min_len_g = $xml_hash->{'qc_params'}{'SRC_min_len'};
    $SRC_min_avg_qual_g = $xml_hash->{'qc_params'}{'SRC_min_avg_qual'};
    $SRC_allow_gaps_g =  $xml_hash->{'qc_params'}{'SRC_allow_gaps'};
    $SRC_allow_ambig_g = $xml_hash->{'qc_params'}{'SRC_allow_ambig'};
    
    #add the samples
    foreach my $sample ( @{$xml_hash->{'sample'}} ) {
        # Guess the name based on the file
        my $name = guess_file_name($sample->{'fwd_file'});
        print "name: $name\n";
        
        # add to listbox
        $list->insert('end', "$name");
        $list->pack();

        push @sample_arr, $sample;
    }

    $mw->update();
}

sub _load_XML_GUI { 
    my ($XML_file_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Load XML File: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 20,
                        -textvariable => $XML_file_ref);
    my $but = $f->Button(-text => "Browse ...",
                         -command => sub { fileDialog($p, $ent, 'open')});
    my $load = $f->Button(-text => "Load",
                           -command => sub { load_XML() });
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left',-expand => 'yes', -fill => 'x');
    $but->pack(-side => 'left');
    $load->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

# method to save current config to be used when the user forgets to add flash/etc ... for the 100th time.
sub save_current_config {
    
    # translate type to seq_type coding
    my $seq_type_code = _translate_seq_type($seq_type_g);
    
    # set up the xml data hash ref    
    my $xml_href = {
        output_root => $output_root_g,
		split_samples_root_dir => '',
        parallelize => $parallelize_bool_g,
        fwd_primer => $fwd_primer_g,
        rev_primer => $rev_primer_g,
        seq_type => $seq_type_code,
        keep_tmp_files => $keep_tmp_files_g,
        fwd_linker => $fwd_linker_g,
        rev_linker => $rev_linker_g,
        fwd_mer_len => $fwd_mer_len_g,
        rev_mer_len => $rev_mer_len_g,
        fwd_max_shifts => $fwd_max_shifts_g,
        rev_max_shifts => $rev_max_shifts_g,
        min_con_depth => $min_con_depth_g,
        diginorm_max => $diginorm_max_g,
        con_algo => $con_algo_g,
        flash_params => {
            m => $m_g,
            M => $M_g,
            x => $x_g,
            p => $p_g,
            r => $r_g,
            f => $fr_g,
            s => $s_g,
        },
        qc_params => {
            con_trim_to_base => $con_trim_to_base_g,
            con_min_len => $con_min_len_g,
            con_min_avg_qual => $con_min_avg_qual_g,
            con_allow_ambig => $con_allow_ambig_g,
            con_min_c_score => $con_min_c_score_g,
            SRC_trim_to_base => $SRC_trim_to_base_g,
            SRC_min_len => $SRC_min_len_g,
            SRC_min_avg_qual => $SRC_min_avg_qual_g,
            SRC_allow_gaps => $SRC_allow_gaps_g,
            SRC_allow_ambig => $SRC_allow_ambig_g,
        },
        sample => \@sample_arr,
    };

    # Build a params object
    my $params_obj = MTToolbox::MTParams->new({href => $xml_href});
    
    # print the config file
    warn("\nSaving temporary config file as './temp_config.xml'. \nWARNING!! This file hasn't been checked for errors.\n\n");
    my $config_file = "temp_config.xml"; 
    $params_obj->print_xml_file($config_file);
}
# </hc>

sub _output_root_field {
    my ($output_root_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Output Directory Root: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 20,
                        -textvariable => $output_root_ref);
    my $but = $f->Button(-text => "Browse ...",
                         -command => sub { dirDialog($p, $ent)});
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left',-expand => 'yes', -fill => 'x');
    $but->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _parallelize_field {
    my ($parallelize_ref, $p) = @_;
        
    my $f = $p->Frame;
    $f->Label(-text => "Parallelize: ")->pack(-side => 'left');
    foreach my $option ( qw/Y N/ ) {
        $f->Radiobutton(-text => $option,
                         -variable => $parallelize_ref,
                         -value => $option
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _fwd_read_file_field {
    my ($fwd_file_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "FWD Reads File: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 20,
                        -textvariable => $fwd_file_ref);
    my $but = $f->Button(-text => "Browse ...",
                         -command => sub { fileDialog($p, $ent, 'open')});
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left',-expand => 'yes', -fill => 'x');
    $but->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _rev_read_file_field {
    my ($rev_file_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "REV Reads File: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 20,
                        -textvariable => $rev_file_ref);
    my $but = $f->Button(-text => "Browse ...",
                         -command => sub { fileDialog($p, $ent, 'open')});
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left',-expand => 'yes', -fill => 'x');
    $but->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _fwd_primer_field {
    my ($fwd_primer_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "FWD Primer: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 50,
                         -textvariable => $fwd_primer_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _rev_primer_field {
    my ($rev_primer_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "REV Primer: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 50,
                         -textvariable => $rev_primer_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _seq_type_field {
    my ($seq_type_ref, $p) = @_;
    
    # 1 => Single End (SE)
    # 2 => Paired End w/ Overlap (Overlapping PE)
    # 3 => Paired End w/o Overlap (Non-Overlapping PE)
    
    my $f = $p->Frame;
    $f->Label(-text => "Read Type: ")->pack(-side => 'left');
    foreach my $option ( ("SE", "PE w/ Overlap", "PE w/o Overlap") ) {
        $f->Radiobutton(-text => $option,
                         -variable => $seq_type_ref,
                         -value => $option
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _keep_tmp_files_field {
    my ($keep_tmp_files_ref, $p) = @_;
        
    my $f = $p->Frame;
    $f->Label(-text => "Keep Temp Files: ")->pack(-side => 'left');
    foreach my $option ( qw/Y N/ ) {
        $f->Radiobutton(-text => $option,
                         -variable => $keep_tmp_files_ref,
                         -value => $option
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _fwd_mer_len_field {
    my ($fwd_mer_len_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "FWD Mer Length: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $fwd_mer_len_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _rev_mer_len_field {
    my ($rev_mer_len_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "REV Mer Len: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $rev_mer_len_ref,
                        )->pack(-side => 'left');

    $f->pack(-fill => 'x');
}

sub _fwd_max_shifts_field {
    my ($fwd_max_shifts_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "FWD Max Frameshifts: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $fwd_max_shifts_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _rev_max_shifts_field {
    my ($rev_max_shifts_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "REV Max Frameshifts: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $rev_max_shifts_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _fwd_linker_field {
    my ($fwd_linker_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "FWD Linker: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 10,
                         -textvariable => $fwd_linker_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _rev_linker_field {
    my ($rev_linker_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "REV Linker: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 10,
                         -textvariable => $rev_linker_ref,
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _min_con_depth_field {
    my ($min_con_depth_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min Con Depth: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $min_con_depth_ref);
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _diginorm_max_field {
    my ($diginorm_max_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Digital Norm Max: ",
                        -anchor => 'e');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $diginorm_max_ref);
    $lab->pack(-side => 'left');
    $ent->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _con_algo_field {
    my ($con_algo_ref, $p) = @_;
    
    my $f = $p->Frame();
    $f->Label(-text => "Consensus Algorithm: ")->pack(-side => 'left');
    foreach my $option ( ("Muscle", "Clustalw", "NoMSA") ) {
        $f->Radiobutton(-text => $option,
                         -variable => $con_algo_ref,
                         -value => $option,
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _min_overlap_field {
    my ($m_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Min Overlap: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $m_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _max_overlap_field {
    my ($M_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Max Overlap: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $M_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _mismatch_ratio_field {
    my ($x_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Mismatch Ratio: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $x_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _phred_offset_field {
    my ($p_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Phred Offset: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $p_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _avg_read_len_field {
    my ($r_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Avg Read Length: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $r_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _avg_frag_len_field {
    my ($fr_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "Avg Fragment Length: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $fr_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _sd_field {
    my ($s_ref, $p) = @_;
    
    my $f = $p->Frame;
    my $lab = $f->Label(-text => "SD of Fragments: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                         -textvariable => $s_ref
                        )->pack(-side => 'left');
    
    $f->pack(-fill => 'x');
}

sub _con_trim_to_base_field {
    my ($con_trim_to_base_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Trim to position: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $con_trim_to_base_ref
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _con_min_len_field {
    my ($con_min_len_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min Length: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $con_min_len_ref
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _con_min_avg_qual_field {
    my ($con_min_avg_qual_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min Average Qual: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $con_min_avg_qual_ref,
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _con_allow_ambig_field {
    my ($con_allow_ambig_ref, $p) = @_;
    
    my $f = $p->Frame;
    $f->Label(-text => "Allow Ambiguous Bases: ")->pack(-side => 'left');
    foreach my $option ( ("Y", "N") ) {
        $f->Radiobutton(-text => $option,
                         -variable => $con_allow_ambig_ref,
                         -value => $option,
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _con_min_c_score_field {
    my ($con_min_c_score_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min c-score: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $con_min_c_score_ref,
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _SRC_trim_to_base_field {
    my ($SRC_trim_to_base_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Trim to position: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $SRC_trim_to_base_ref
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _SRC_min_len_field {
    my ($SRC_min_len_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min Length: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $SRC_min_len_ref
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _SRC_min_avg_qual_field {
    my ($SRC_min_avg_qual_ref, $p) = @_;
    
    my $f = $p->Frame();
    my $lab = $f->Label(-text => "Min Average Qual: ")->pack(-side => 'left');
    my $ent = $f->Entry(-width => 5,
                        -textvariable => $SRC_min_avg_qual_ref,
                        )->pack(-side => 'left');
    $f->pack(-fill => 'x');
}

sub _SRC_allow_gaps_field {
    my ($SRC_allow_gaps_ref, $p) = @_;
    
    my $f = $p->Frame;
    $f->Label(-text => "Allow Gaps: ")->pack(-side => 'left');
    foreach my $option ( ("Y", "N") ) {
        $f->Radiobutton(-text => $option,
                         -variable => $SRC_allow_gaps_ref,
                         -value => $option,
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

sub _SRC_allow_ambig_field {
    my ($SRC_allow_ambig_ref, $p) = @_;
    
    my $f = $p->Frame;
    $f->Label(-text => "Allow Ambiguous Bases: ")->pack(-side => 'left');
    foreach my $option ( ("Y", "N") ) {
        $f->Radiobutton(-text => $option,
                         -variable => $SRC_allow_ambig_ref,
                         -value => $option,
                         )->pack(-side => 'left');
    }
    $f->pack(-fill => 'x');
}

__END__

# POD

=head1 NAME

MTToolbox.pl - Runs the GUI for MTToolbox analysis


=head1 VERSION

This documentation refers to MTToolbox.pl version 4.1.2


=head1 USAGE

MTToolbox.pl


=head1 REQUIRED ARGUMENTS

NA
    
    
=head1 OPTIONAL ARGUMENTS

NA


=head1 DESCRIPTION

Limitting the number of sequencing errors is crutial for precise sequence
comparisons.  Sequencing errors are often introduced during the PCR steps
and/or the sequencing step.  In an attempt to limit the number of errors
we have developed a molecular biology technique in which individual DNA
molecules are tagged with a unique barcode-like sequences.  We call this
sequence a molecule tag (MT).  After each DNA molecule has been tagged it
can be PCRed and sequenced.  After sequencing we recieve raw reads that
contain a molecule tag at the beginning.  MTToolbox.pl sorts those sequences
based on their molecule tags and builds a consensus sequence to represent
the original DNA molecule that was tagged.

For large sequence sets across several samples, these operations can be quite
time intensive so MTToolbox.pl includes an option to run in parallel using LSF.  


=head1 OUTPUTS

The following is a list and summary of some of the major output files:

all_summary.txt - Summary numbers for the entire run and each sample
individually.

all_summary.ps - PostScript version of the above text file.

config.xml - The config file with all the run variables.

all_consensusSeqs_LQ.fastq - Consensus sequences that failed the QC test.

all_consensusSeqs_diagnostics.txt - Reasons why each of the low quality reads
failed the QC tests.  If large numbers of sequences failed the QC tests this
file can be used to diagnose what went wrong.

all_single_read_categories_HQ.fastq - MT categories with only one raw read.  It
is imposible to build a consensus sequences from just one read, however it these
may still be useful in downstream analysis.

all_single_read_categories_LQ.fastq - MT categories with only one raw read that
were classified as low quality

all_single_read_categories_diagnostics.txt - Reasons why each of the low quality
reads failed the QC tests.

all_categorizable_reads.fastq - All raw reads which match the expected pattern.

all_seqs_per_mt_dist.txt - Distribution of raw reads per MT category text file

all_seqs_per_mt_hist.txt - The above text file in histogram format

all_seqs_per_mt_hist.ps - The PostScript version of the above file

con_ambig_base_comp_dist.png - Graph of the kinds of ambiguous bases in the
consensus seqs

con_ambig_base_pos_dist.png - Graph of where the ambiguous base position in the
consensus seqs.

con_ambig_base_per_seq.png - Graph showing the distribution of ambiguous bases
in each read.
    

=head1 CONFIGURATION AND ENVIRONMENT
    
Tested using perl5.8.9 and perl5.12.3.
    
    
=head1 DEPENDANCIES

    Tk
    Tk::NoteBook
    Tk::DialogBox
    Tk::LabEntry
    XML::Simple
    Carp
    Data::Dumper
    MTToolbox::MTParams 4.1.2
    MTToolbox::MTEnv 4.1.2
    
    Note: Perl Tk requires a XWindow server
    

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to author


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


=cut
